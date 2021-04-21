#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <algorithm>

using namespace std;

struct Object {
	int cost = 0;
	int weight = 0;
	bool used = false;
	Object(int obj_cost, int obj_weight) {
		cost = obj_cost;
		weight = obj_weight;
	}
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

int main() {
	vector<Object> objects = {
		Object(4, 3),
		Object(2, 3),
		Object(1, 2),
		Object(3, 5),
		Object(5, 6),
		Object(2, 2),
		Object(9, 8),
	};

	unsigned int start_time = clock(); // начальное время
	int cost = MDP2(objects, 20);
	unsigned int end_time = clock(); // конечное время
	
	cout << "cost:" << cost << endl;
	cout << "time:" << end_time - start_time << endl; // искомое время
	return 0;



	cout << "C best " << ptas(objects, 11, 3) << endl;
	return 0;
}