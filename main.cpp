#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <algorithm>
#include <fstream>

using namespace std;

struct Object {
    int id;
	int cost = 0;
	int weight = 0;
	double quality = 0.0;
	bool used = false;
	short use_condition = -1;
	Object(int obj_cost, int obj_weight, int obj_id) :
            cost(obj_cost), weight(obj_weight), id(obj_id){
	};
};

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

int greedyAlgorithm(vector<Object>& objects, int capacity){
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

int maxGreed(vector<Object>& objects, int capacity){
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

vector<float> relaxedGreedyAlgorithm(vector<Object>& objects, double capacity){
    vector<float> res(objects.size());
    for(auto& obj : objects) {
        if(obj.use_condition == 1) {
            capacity -= obj.weight;
            res[obj.id] = 1;
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
            } else {
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

int main() {

    vector<Object> objj = {
            {6, 1, 0},
            {10, 2, 1},
            {10, 3, 2}
    };
    branchAndBound(objj, 5);
    auto ttmp = totalcost;
    for(auto n : totalvec)
        cout<<n<<" ";
    cout<<endl;

    string file_template="../data/";

    int capacity = 0;
    int p = 0;
    int w = 0;
    vector<Object> objects;

    for(int i=1; i<8; ++i) {
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

        auto res = twoApporxAlgorithm(objects, capacity);
        for(auto n : res)
            cout<<n<<" ";
        cout<<endl;
        totalvec.clear();
        totalcost = 0;
     //   branchAndBound(objects, capacity);
        cout<<totalcost<<endl;
        for(auto n : totalvec)
            cout<<n<<" ";
        cout<<endl;
        // TO DO INSERT ALGOS HERE
        objects.clear();
    }

    return 0;
}