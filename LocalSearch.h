#pragma once

#include "SymmetricMatrix.h"

#include <random>
#include <string>

using namespace std;

class QAP_data
{
public:
	QAP_data(string filename);
	void out_dists(ostream& out);
	void out_flows(ostream& out);
	long long  fitness_function(vector<int>& location2factory);
	int size();

private:
	SymmetricMatrix dists;
	SymmetricMatrix flows;
};

QAP_data::QAP_data(string filename) {
	ifstream fin;
	fin.open("../data/tai20a");
	if (!fin) {
		throw - 1;
	}

	int size;
	fin >> size;

	dists.read_from_taiNa(fin, size);
	flows.read_from_taiNa(fin, size);

	fin.close();
}

void QAP_data::out_dists(ostream& out) {
	out << dists << endl;
}

void QAP_data::out_flows(ostream& out) {
	out << flows << endl;
}

long long  QAP_data::fitness_function(vector<int>& location2factory) {
	int size = dists.size();
	long long fitness_value = 0;

	for (int location_i = 0; location_i < size; ++location_i) {
		for (int location_j = location_i; location_j < size; ++location_j) {
			long long dist = dists.item(location_i, location_j);
			int factory_i = location2factory[location_i];
			int factory_j = location2factory[location_j];
			long long flow = flows.item(factory_i, factory_j);
			fitness_value += dist * flow;
		}
	}
	return fitness_value;
}

int QAP_data::size() {
	return dists.size();
}



vector<int>& BestLocalSearch(vector<int>& location2factory, QAP_data& data, 
	long long& fitness_function_value) {
	vector<int> location2factory_buf = location2factory;
	long long ffunc_value_best = data.fitness_function(location2factory);
	int size = location2factory.size();

	for (int i = 0; i < size - 1; ++i) {
		for (int j = i + 1; j < size; ++j) {
			// stochastic 2-opt
			int first = rand() % size;
			int second = rand() % size;
			if (second < first) swap(first, second);

			int first_copy = first;
			int second_copy = second;

			while (first_copy < second_copy) {
				swap(location2factory_buf[first_copy], location2factory_buf[second_copy]);
				++first_copy;
				--second_copy;
			}

			//swap(location2factory_buf[i], location2factory_buf[j]);

			long long new_ffunc_value = data.fitness_function(location2factory_buf);
			if (new_ffunc_value < ffunc_value_best) {
				location2factory = location2factory_buf;
				ffunc_value_best = new_ffunc_value;
			}

			//swap(location2factory_buf[i], location2factory_buf[j]);
			while (first < second) {
				swap(location2factory_buf[first], location2factory_buf[second]);
				++first;
				--second;
			}
			int a = 5;
		}
	}

	fitness_function_value = ffunc_value_best;
	return location2factory;
}

void perturbation(vector<int>& vec, int pert_coef) {
	int size = vec.size();
	for (int i = 0; i < pert_coef; ++i) {
		int first = rand() % size;
		int second = rand() % size;

		swap(vec[first], vec[second]);
	}
}

vector<int> IteratedLocalSearch(string tai_filename, int iterations, int pert_coef,
	long long& best_value) {
	// getting the data
	QAP_data data(tai_filename);
	//data.out_dists(cout);
	//data.out_flows(cout);

	// generating a simple solution
	int size = data.size();
	vector<int> location2factory(size);
	for (int i = 0; i < size; ++i) {
		location2factory[i] = i;
	}

	// local search
	long long fitness_function_value = 0;  // just an additional returned value
	location2factory = BestLocalSearch(location2factory, data, fitness_function_value);
	vector<int> location2factory_best = location2factory;

	for (int i = 0; i < iterations; ++i) {
		// perturbation
		perturbation(location2factory, pert_coef);

		// local search
		long long fitness_function_value_perturb = 0;
		vector<int>& location2factory_perturb = BestLocalSearch(location2factory, data, 
			fitness_function_value_perturb);

		// accept best
		if (fitness_function_value_perturb < fitness_function_value) {
			location2factory_best = location2factory_perturb;
			location2factory = location2factory_perturb;
			fitness_function_value = fitness_function_value_perturb;
		}
	}

	best_value = fitness_function_value;
	return location2factory;
}
