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

int ComputeCost(vector<int> ids, vector<Object>& objects, int W) {
	int current_W = 0;
	int current_C = 0;
	for (int& id : ids) {
		current_W += objects[id].weight;
		current_C += objects[id].cost;
		objects[id].used = true;
	}
	if (current_W > W) {
		return 0;
	}
	for (int i = 0; i < objects.size(); i++) {

		if (objects[i].used) continue;
		int new_W = current_W + objects[i].weight;
		int new_C = current_C + objects[i].cost;
		if (new_W > W) {
			return current_C;
		}
		else {
			current_W = new_W;
			current_C = new_C;
		}
	}
	for (auto& object : objects) {
		object.used = false;
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

int ComputeCBest(vector<Object>& objects, int W, int k) {
	vector<int> ids(k);
	int value = 0;
	for (int& id : ids) {
		id = value;
		value++;
	}

	int C_best = 0;
	if (k == 0) {
		return ComputeCost(ids, objects, W);
	}
	
	while (true) {
		int C_new = ComputeCost(ids, objects, W);
		C_best = C_new > C_best ? C_new : C_best;

		int obj_number = objects.size();

		if (RecountIds(ids, obj_number))
			return C_best;
	}
}

int ptas(vector<Object>& objects, int W, int k) {
	int C_best = 0;
	
	for (int i = 0; i <= k; i++) {
		int C_new = ComputeCBest(objects, W, i);
		C_best = C_new > C_best ? C_new : C_best;
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
	vector<int>& prev, vector<int>& curr, const vector<Object>& objects, int& W) {

	for (int item_index = start_index; item_index < finish_index; ++item_index) {
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
}

int MDP2_speedup_nt(vector<Object>& objects, int W) {
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

	for (int string_index = 1; string_index <= size; ++string_index) {
		start_index = max(barrier[string_index], 0);
		if (false) {
			mid = (start_index + length) / 2;
			thread left_work(compute_curr, ref(start_index), ref(mid), ref(string_index),
				ref(prev), ref(curr), ref(objects), ref(W));
			compute_curr(mid, length, string_index, prev, curr, objects, W);
			left_work.join();
		}
		else {
			//cout << start_index << ":" << length << endl;
			compute_curr(start_index, length, string_index, prev, curr, objects, W);
		}
		swap(prev, curr);
	}
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
			compute_curr(mid, length, string_index, prev, curr, objects, W);

			std::unique_lock<std::mutex> ul(m);

			cv.wait(ul, [=]() { return NewComp == false; });

			ul.unlock();
		}
		else {
			//cout << start_index << ":" << length << endl;
			compute_curr(start_index, length, string_index, prev, curr, objects, W);
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
    sort(objects.begin(), objects.end(), [](Object& objl, Object& objr){
        return objl.quality > objr.quality;
    });
    for(size_t i = 0; i<objects.size(); ++i){
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
    sort(objects.begin(), objects.end(), [](Object& objl, Object& objr){
        return objl.cost > objr.cost;
    });
    for(size_t i = 0; i<objects.size(); ++i){
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
    if(greedyAlgorithm(tmp1, capacity) > maxGreed(tmp2, capacity)){
        for(auto& obj: tmp1)
            if(obj.used) res[obj.id] = true;
    }
    else{
        for(auto& obj: tmp2){
            if(obj.used) res[obj.id] = true;
        }
    }
    return res;
}

vector<float> relaxedGreedyAlgorithm(vector<Object> objects, double capacity){
    vector<float> res(objects.size(), 0);
    for(auto& obj : objects) {
        if(obj.use_condition == 1) {
            if (capacity - obj.weight >= 0) {
                capacity -= obj.weight;
                res[obj.id] = 1;
                }
            else return res;
        }
        else
            obj.quality = ((double) obj.cost / (double) obj.weight);
    }
    sort(objects.begin(), objects.end(), [](Object& objl, Object& objr){
        return objl.quality > objr.quality;
    });
    for(size_t i = 0; i<objects.size(); ++i) {
        if (objects[i].use_condition==-1) {
            if (capacity - objects[i].weight >= 0) {
                capacity -= objects[i].weight;
                res[objects[i].id] = 1;
            } else if(capacity>=0) {
                float a = (float) capacity / (float) objects[i].weight;
                capacity -= a * (float)objects[i].weight;
                res[objects[i].id] = a;
            }
        }
    }
    return res;
}

vector<float> totalvec;
int totalcost=0;

void branchAndBound(vector<Object> objects, int capacity){
    vector<float> resvec = relaxedGreedyAlgorithm(objects, capacity);
    int tmpcost = 0;
    for(size_t i = 0; i<objects.size(); ++i)
        if(resvec[i]) tmpcost+=objects[i].cost*resvec[i];
    if(tmpcost>totalcost) {
        for (size_t i = 0; i < objects.size(); ++i)
            if (resvec[i] < 1 && resvec[i] > 0) {
                objects[i].use_condition = 1;
                branchAndBound(objects, capacity);
                objects[i].use_condition = 0;
                branchAndBound(objects, capacity);
                return;
            }
    }
    if(tmpcost>totalcost){
        totalvec=resvec;
        totalcost=tmpcost;
    }
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
	std::cout << "Sort-barrier + 2-proc speedup (only for really large cases): " << cs << "$. The time: " << elapseds_ms.count() << " ms\n";
}

int main() {
	//MDP2_comp();
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
                  << objects.size()*log(objects.size())<<"\n"
                  << std::endl;

        totalvec.clear();
        totalcost = 0;
        cout<<"branchAndBound answer: \n";
        branchAndBound(objects, capacity);
        for(auto n : totalvec)
            cout<<n<<" ";
        tmp = get_cost_and_weight<int>(totalvec, objects);
        cout<<", total cost = "<<get<0>(tmp)<<" total weight "<<get<1>(tmp)<<endl;
        pr_StartTime = std::chrono::steady_clock::now();
        global_counter = 0;
        for(size_t i = 0; i < 10; ++i) {
            totalvec.clear();
            totalcost = 0;
            branchAndBound(objects, capacity);
        }
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10
                  << " number of comparisons (used estimated number of comparisons in quick sort (nlogn)) = "
                  << objects.size()*log(objects.size())<<"\n"
                  << std::endl;
        // TO DO INSERT ALGOS HERE
        objects.clear();
    }

    return 0;
}