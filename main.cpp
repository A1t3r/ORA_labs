#include "LocalSearch.h"

#include <iostream>

using namespace std;

int main() {
	vector<int> location2factory;
	long long best_value = -1;

	location2factory = IteratedLocalSearch("../data/tai20a", 20, best_value);
	cout << "Value " << best_value << endl << endl;
	for (int i = 0; i < location2factory.size(); ++i) {
		cout << location2factory[i] << " " << i << endl;
	}
	cout << endl;

	return 0;
}