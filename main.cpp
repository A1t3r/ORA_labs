#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <utility>
#include <thread>
#include <mutex>

#include "Object.h"

using namespace std;

int global_counter;

mutex m;
condition_variable cv;
bool NewComp = false;
bool FinishThread = false;

template<class T>
pair<T, T> get_cost_and_weight(vector<bool>& res, vector<Object>& objs){
    T tmpcost = 0;
    T tmpweight = 0;
    for(size_t i = 0; i< objs.size(); ++i){
        if(res[i]){
            tmpcost += objs[i].cost;
            tmpweight += objs[i].weight;
        }
    }
    return{tmpcost, tmpweight};
}

template<class T>
pair<T, T> get_cost_and_weight(vector<float>& res, vector<Object>& objs){
    T tmpcost = 0;
    T tmpweight = 0;
    for(size_t i = 0; i< objs.size(); ++i){
        if(res[i]){
            tmpcost += objs[i].cost;
            tmpweight += objs[i].weight;
        }
    }
    return{tmpcost, tmpweight};
}

int ptas_W = 0;
vector<int> ptas_ids = {};

int ComputeCost(vector<int> ids, vector<Object>& objects, int W, int& buf_w, vector<int>& buf_ids) {
	int current_W = 0;
	int current_C = 0;
	for (int& id : ids) {
		current_W += objects[id].weight;
		current_C += objects[id].cost;
		objects[id].used = true;
	}

	buf_w = current_W;
	buf_ids = ids;

	if (current_W > W) {
		buf_w = 0;
		buf_ids = {};
		return 0;
	}
	for (int i = 0; i < objects.size(); i++) {

		if (objects[i].used) continue;
		int new_W = current_W + objects[i].weight;
		int new_C = current_C + objects[i].cost;
		objects[i].used = true;
		if (new_W > W) {
			return current_C;
		}
		else {
			current_W = new_W;
			current_C = new_C;

			buf_w = new_W;
		}
	}

	buf_ids = {};

	for (int i = 0; i < objects.size(); i++) {
		if (objects[i].used) {
			objects[i].used = false;
			buf_ids.push_back(i);
		}
	}

	return current_C;
}

bool NextSet(vector<int>& a, int n, int m) {
	int k = m;
	for (int i = k - 1; i >= 0; --i)
		if (a[i] < n - k + i + 1)
		{
			++a[i];
			for (int j = i + 1; j < k; ++j)
				a[j] = a[j - 1] + 1;
			return true;
		}
	return false;
}

bool RecountIds(vector<int>& ids, int obj_number) {
	return !NextSet(ids, obj_number - 1, ids.size());
}

int ComputeCBest(vector<Object>& objects, int W, int k, int& buf_w, vector<int>& buf_ids) {
	vector<int> ids(k);
	int value = 0;
	for (int& id : ids) {
		id = value;
		value++;
	}

	int ccbbuf_w = 0;
	vector<int> ccbbuf_ids = {};

	int C_best = 0;
	if (k == 0) {
		return ComputeCost(ids, objects, W, buf_w, buf_ids);
	}

	while (true) {
		int C_new = ComputeCost(ids, objects, W, ccbbuf_w, ccbbuf_ids);

		C_best = C_new > C_best ? C_new : C_best;
		if (C_new > C_best) {
			C_best = C_new;
			buf_w = ccbbuf_w;
			buf_ids = ccbbuf_ids;
		}

		int obj_number = objects.size();

		if (RecountIds(ids, obj_number))
			return C_best;
	}
}

int ptas(vector<Object>& objects, int W, int k) {
	int C_best = 0;
	int buf_w = 0;
	vector<int> buf_ids = {};
	ptas_ids = {};
	ptas_W = 0;
	
	for (int i = 0; i <= k; i++) {
		int C_new = ComputeCBest(objects, W, i, buf_w, buf_ids);
		if (C_new > C_best) {
			C_best = C_new;
			ptas_ids = buf_ids;
			ptas_W = buf_w;
		}
	}
	return C_best;
}

int MDP2(vector<Object>& objects, int W) {
	int length = W + 1;
	int size = objects.size();
	vector<int> prev(length);
	vector<int> curr(length);
	vector<int> barrier(size + 1);
	//barrier[size] = length - 1;
	// step 1: computing barrier
	//for (int i = size - 1; i >= 0; --i){
	//	barrier[i] = barrier[i + 1] - objects[i].weight;
	//}

	for (int string_index = 1; string_index <= size; ++string_index) {
		for (int item_index = max(barrier[string_index], 0); item_index < length; ++item_index) {
			int obj_W = objects[string_index - 1].weight;
			if (obj_W <= W) {
				int prev_bp = 0;
				if (item_index >= obj_W) {
					prev_bp = objects[string_index - 1].cost + prev[item_index - obj_W];
				}

				curr[item_index] = max(
					prev[item_index], 
					prev_bp
					);
			}
			else {
				break;
			}
		}
		swap(prev, curr);
	}
	return prev[length - 1];
}

void compute_curr_thread(int& start_index, int& finish_index, int& string_index,
	vector<int>& prev, vector<int>& curr, const vector<Object>& objects, int& W) {
	vector<int>& thread_curr = curr;
	vector<int>& thread_prev = prev;

	while (true) {
		std::unique_lock<mutex> ul(m);
		cv.wait(ul, [=]() { return NewComp; });

		if (FinishThread) break;


		//cout << "thread " << start_index << ":" << finish_index << endl;
		for (int item_index = start_index; item_index < finish_index; ++item_index) {
			int obj_W = objects[string_index - 1].weight;
			if (obj_W <= W) {
				int prev_bp = 0;
				if (item_index >= obj_W) {
					prev_bp = objects[string_index - 1].cost + thread_prev[item_index - obj_W];
				}

				thread_curr[item_index] = max(
					thread_prev[item_index],
					prev_bp
				);
			}
			else {
				break;
			}
		}

		NewComp = false;
		cv.notify_one();
	}
}

void compute_curr(int& start_index, int& finish_index, int& string_index,
	vector<int>& prev, vector<int>& curr, const vector<Object>& objects, int& W,
	vector<vector<int>>& history) {

	vector<vector<int>> new_history(history.size());

	for (int item_index = start_index; item_index < finish_index; ++item_index) {
		int obj_W = objects[string_index - 1].weight;
		if (obj_W <= W) {
			int prev_bp = 0;
			if (item_index >= obj_W) {
				prev_bp = objects[string_index - 1].cost + prev[item_index - obj_W];
			}

			if (prev[item_index] >= prev_bp) {
				curr[item_index] = prev[item_index];
				new_history[item_index] = history[item_index];
			}
			else {
				curr[item_index] = prev_bp;
				new_history[item_index] = history[item_index - obj_W];
				new_history[item_index].push_back(string_index - 1);
			}
		}
		else {
			break;
		}
	}
	swap(new_history, history);
}

vector<int> MPD_taken = {};

int MDP2_speedup_nt(vector<Object>& objects, int W) {
	if (objects.size() == 0) {
		return 0;
	}

	sort(objects.begin(), objects.end(), [](const Object& lhs, const Object& rhs) -> bool
		{
			return lhs.weight > rhs.weight;
		});
	int length = W + 1;
	int size = objects.size();
	int mid = 0;
	int start_index = 0;

	vector<int> prev(length);
	vector<int> curr(length);
	vector<int> barrier(size + 1);
	int threshold = 100;
	barrier[size] = length - 1;
	// step 1: computing barrier
	for (int i = size - 1; i >= 0; --i) {
		barrier[i] = barrier[i + 1] - objects[i].weight;
	}

	vector<vector<int>> history(length);

	for (int string_index = 1; string_index <= size; ++string_index) {
		start_index = max(barrier[string_index], 0);
		compute_curr(start_index, length, string_index, prev, curr, objects, W, history);
		swap(prev, curr);
	}

	MPD_taken = history[length - 1];
	return prev[length - 1];
}

int MDP2_speedup(vector<Object>& objects, int W) {
	sort(objects.begin(), objects.end(), [](const Object& lhs, const Object& rhs) -> bool
		{
			return lhs.weight > rhs.weight;
		});
	int length = W + 1;
	int size = objects.size();
	int mid = 0;
	int start_index = 0;
	int obj_id = 0;

	vector<int> prev(length);
	vector<int> curr(length);
	vector<int> barrier(size + 1);
	int threshold = 100;
	barrier[size] = length - 1;
	// step 1: computing barrier
	for (int i = size - 1; i >= 0; --i){
		barrier[i] = barrier[i + 1] - objects[i].weight;
	}

	vector<vector<int>> history = {};

	thread left_work(compute_curr_thread, ref(start_index), ref(mid), ref(obj_id),
		ref(prev), ref(curr), ref(objects), ref(W));

	for (int string_index = 1; string_index <= size; ++string_index) {
		//cout << string_index << endl;
		obj_id = string_index;
		start_index = max(barrier[string_index], 0);
		if (length - start_index > threshold) {
			mid = (start_index + length) / 2;
			{
				std::unique_lock<std::mutex> ul(m);
				NewComp = true;
				cv.notify_one();
			}
			//cout << "main " << mid << ":" << length << endl;
			compute_curr(mid, length, string_index, prev, curr, objects, W, history);

			std::unique_lock<std::mutex> ul(m);

			cv.wait(ul, [=]() { return NewComp == false; });

			ul.unlock();
		}
		else {
			//cout << start_index << ":" << length << endl;
			compute_curr(start_index, length, string_index, prev, curr, objects, W, history);
		}
		swap(prev, curr);
	}
	FinishThread = true;
	NewComp = true;
	cv.notify_all();
	left_work.join();

	return prev[length - 1];
}

static int greedyAlgorithm(vector<Object>& objects, int capacity){
    int res = 0;
    for(auto& obj : objects)
        obj.quality = ((double)obj.cost / (double)obj.weight);
    global_counter++;
    sort(objects.begin(), objects.end(), [](Object& objl, Object& objr){
        return objl.quality > objr.quality;
    });
    for(size_t i = 0; i<objects.size(); ++i){
        global_counter++;
        if(capacity - objects[i].weight >=0){
            capacity -= objects[i].weight;
            objects[i].used = true;
            res += objects[i].cost;
        }
    }
    return res;
}

static int maxGreed(vector<Object>& objects, int capacity){
    int res = 0;
    global_counter++;
    sort(objects.begin(), objects.end(), [](Object& objl, Object& objr){
        return objl.cost > objr.cost;
    });
    for(size_t i = 0; i<objects.size(); ++i){
        global_counter++;
        if(capacity - objects[i].weight >=0){
            capacity -= objects[i].weight;
            objects[i].used = true;
            res += objects[i].cost;
        }
    }
    return res;
}

vector<bool> twoApporxAlgorithm(vector<Object>& objects, int capacity){
    vector<bool> res(objects.size(), false);
    auto tmp1 = objects;
    auto tmp2 = objects;
    global_counter++;
    if(greedyAlgorithm(tmp1, capacity) > maxGreed(tmp2, capacity)){
        for(auto& obj: tmp1)
            if(obj.used) res[obj.id] = true;
    }
    else{
        for(auto& obj: tmp2){
            if(obj.used) res[obj.id] = true;
        }
    }
    global_counter+=tmp1.size();
    return res;
}

double relaxedGreedyAlgorithm(vector<Object>& objects, double capacity){
    double tmpcost = 0;
    for(auto& obj : objects) {
     //   global_counter++;
        if(obj.use_condition == 1) {
         //   global_counter++;
            if (capacity - obj.weight >= 0) {
                tmpcost += obj.cost;
                capacity -= obj.weight;
                obj.used_scale = 1;
                }
            else return tmpcost;
        }
    }
    for(size_t i = 0; i<objects.size(); ++i) {
        if (objects[i].use_condition==-1) {
       //     global_counter++;
            if (capacity - objects[i].weight >= 0) {
                capacity -= objects[i].weight;
                objects[i].used_scale = 1;
                tmpcost += objects[i].cost;
            } else if(capacity>=0) {
          //      global_counter++;
                float a = (float) capacity / (float) objects[i].weight;
                capacity -= a * (float)objects[i].weight;
                objects[i].used_scale = a;
                tmpcost += objects[i].cost * a;
            }
        }
    }
    return tmpcost;
}


void branchAndBound(vector<Object> objects, int capacity, double& totalcost, vector<float>& totalvec) {
    global_counter++;
    auto tmpobjects = objects;
    double tmpcost = relaxedGreedyAlgorithm(tmpobjects, capacity);
    if (tmpcost > totalcost) {
        for (size_t i = 0; i < objects.size(); ++i) {
          //  global_counter++;
            if (tmpobjects[i].used_scale < 1 && tmpobjects[i].used_scale > 0) {
                objects[i].use_condition = 1;
                branchAndBound(objects, capacity, totalcost, totalvec);
                objects[i].use_condition = 0;
                branchAndBound(objects, capacity, totalcost, totalvec);
                return;
            }
        }
        if (tmpcost > totalcost) {
          //  global_counter++;
            for (auto item : tmpobjects)
                if (item.used_scale){
          //          global_counter++;
                    totalvec[item.id] = 1;
                }
                else totalvec[item.id] = 0;
            totalcost = tmpcost;
        }
    }
}

vector<float> BBM(vector<Object> objects, int capacity){
    vector<float> totalvec(objects.size(), 0);
    double totalcost=0;
    for(auto& obj : objects)
        obj.quality = ((double) obj.cost / (double) obj.weight);
    sort(objects.begin(), objects.end(), [](Object& objl, Object& objr){
        return objl.quality > objr.quality;
    });
    branchAndBound(objects, capacity, totalcost,totalvec);
    return totalvec;
}

void MDP2_comp() {
	int W = 18000;

	auto begin = std::chrono::steady_clock::now();
	int c = MDP2(example, W);
	auto end = std::chrono::steady_clock::now();

	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cout << "Base version: " << c << "$. The time: " << elapsed_ms.count() << " ms\n";

	auto beginsnt = std::chrono::steady_clock::now();
	int csnt = MDP2_speedup_nt(example, W);
	auto endsnt = std::chrono::steady_clock::now();

	auto elapsedsnt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endsnt - beginsnt);
	std::cout << "Sort-barrier speedup: " << csnt << "$. The time: " << elapsedsnt_ms.count() << " ms\n";

	auto begins = std::chrono::steady_clock::now();
	int cs = MDP2_speedup(example, W);
	auto ends = std::chrono::steady_clock::now();

	auto elapseds_ms = std::chrono::duration_cast<std::chrono::milliseconds>(ends - begins);
	std::cout << "Sort-barrier + 2-proc speedup: " << cs << "$. The time: " << elapseds_ms.count() << " ms\n";
}

int main() {
	cout << "cost " << ptas(example, 5000, 2) << endl;
	cout << "weight " << ptas_W << endl;
	cout << "Ids:" << endl;
	for (const auto& item : ptas_ids) {
		cout << item << " ";
	}
	cout << endl;
	cout << endl;
	cout << endl;

	cout << "cost " << MDP2_speedup_nt(example, 10) << endl;
	cout << "ids:" << endl;
	int wmpd = 0;
	for (const auto& item : MPD_taken) {
		cout << item << " ";
		wmpd += example[item].weight;
	}
	cout << endl;
	cout << "weight " << wmpd << endl;
	return 0;

    string file_template="../data/";
    std::chrono::steady_clock::time_point pr_StartTime;
    std::chrono::steady_clock::time_point pr_EndTime;
    int capacity = 0;
    int p = 0;
    int w = 0;
    vector<Object> objects;

    for(int i=1; i<8; ++i) {
        cout<<"Now testing sample number "<<to_string(i)<<endl;
        int id = 0;
        ifstream file_capacity(file_template + "p0" + std::to_string(i) + "_c.txt");
        file_capacity>>capacity;
        ifstream file_p(file_template + "p0" + std::to_string(i) + "_p.txt");
        ifstream file_w(file_template + "p0" + std::to_string(i) + "_w.txt");
        while(!file_p.eof()){
            file_p >> p, file_w >> w;
            objects.emplace_back(p, w, id);
            id++;
        }

        cout<<"twoApporxAlgorithm answer: \n";
        vector<bool> res = twoApporxAlgorithm(objects, capacity);
        for(auto n : res)
            cout<<n<<" ";
        auto tmp = get_cost_and_weight<int>(res, objects);
        cout<<", total cost = "<<get<0>(tmp)<<" total weight "<<get<1>(tmp)<<endl;
        pr_StartTime = std::chrono::steady_clock::now();
        global_counter = 0;
        for(size_t i = 0; i < 10; ++i)
            twoApporxAlgorithm(objects, capacity);
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10
                  << " number of comparisons (used estimated number of comparisons in quick sort (nlogn)) = "
                  << (objects.size()*log(objects.size()) * 2 + global_counter) / 10 <<"\n"
                  << std::endl;


        cout<<"branchAndBound answer: \n";
        auto res2 = BBM(objects, capacity);
        for(auto n : res2)
            cout<<n<<" ";
        tmp = get_cost_and_weight<int>(res2, objects);
        cout<<", total cost = "<<get<0>(tmp)<<" total weight "<<get<1>(tmp)<<endl;
        pr_StartTime = std::chrono::steady_clock::now();
        global_counter = 0;
        for(size_t i = 0; i < 10; ++i) {
            BBM(objects, capacity);
        }
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10
                  << " number of branches = "
                  <<  global_counter  / 10 <<"\n"
                  << std::endl;
        // TO DO INSERT ALGOS HERE
        objects.clear();
    }

    return 0;
}