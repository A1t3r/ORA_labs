#pragma once

#include "SymmetricMatrix.h"

#include <random>
#include <string>
#include <omp.h>
#include <iostream>

using namespace std;

class QAP_data
{
public:
	QAP_data(string filename);
	void out_dists(ostream& out);
	void out_flows(ostream& out);
	long long fitness_function(vector<int>& location2factory);
	long long parallel_fitness_function(vector<int>& location2factory);
    long long recalculateFitness(vector<int>& sol, size_t i,
                                           size_t j, long long old_val);
	int size();

private:
	SymmetricMatrix dists;
	SymmetricMatrix flows;
	vector<int> buf;
};

QAP_data::QAP_data(string filename) {
	ifstream fin;
	//fin.open("../data/tai20a");
    fin.open(filename);
	if (!fin) {
		throw - 1;
	}

	int size;
	fin >> size;

	buf.resize(size * size);

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

long long  QAP_data::parallel_fitness_function(vector<int>& location2factory) {
	int size = dists.size();
	long long fitness_value = 0;
	
	
	for (int location_i = 0; location_i < size; ++location_i) {
		for (int location_j = location_i; location_j < size; ++location_j) {
			
			long long dist = dists.item(location_i, location_j);
			int factory_i = location2factory[location_i];
			int factory_j = location2factory[location_j];
			long long flow = dist * flows.item(factory_i, factory_j);

			buf[location_i * size + location_j] = flow;
		}
	}

	for (auto& item : buf) {
		fitness_value += item;
		item = 0;
	}
	return fitness_value;
}

int QAP_data::size() {
	return dists.size();
}

long long QAP_data::recalculateFitness(vector<int>& sol, size_t i, size_t j, long long old_val){
    long long delta = 0;
    for(size_t n = 0; n < sol.size(); ++n){
        delta -= dists.item(i, n) * flows.item(sol[i], sol[n]);
        delta -= dists.item(j, n) * flows.item(sol[j], sol[n]);
    }
    return delta;
}

vector<int>& BestLocalSearch(vector<int>& location2factory, QAP_data& data, 
	long long& fitness_function_value, int opt_coef = 8) {
	vector<int> location2factory_buf = location2factory;
	long long ffunc_value_best = data.fitness_function(location2factory);
	int size = location2factory.size();

	for (int i = 0; i < opt_coef; ++i) {
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

		long long new_ffunc_value = data.fitness_function(location2factory_buf);
		if (new_ffunc_value < ffunc_value_best) {
			location2factory = location2factory_buf;
			ffunc_value_best = new_ffunc_value;
		}

		while (first < second) {
			swap(location2factory_buf[first], location2factory_buf[second]);
			++first;
			--second;
		}
	}

	fitness_function_value = ffunc_value_best;
	return location2factory;
}

vector<int>& BestLocalSearchParallel(vector<int>& location2factory, QAP_data& data,
	long long& fitness_function_value, int opt_coef = 8) {
	long long ffunc_value_best = data.parallel_fitness_function(location2factory);
	int size = location2factory.size();

	for (int i = 0; i < opt_coef; ++i) {
		vector<int> location2factory_buf = location2factory;
		
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

		long long new_ffunc_value = data.parallel_fitness_function(location2factory_buf);
		if (new_ffunc_value < ffunc_value_best) {
			location2factory = location2factory_buf;
			ffunc_value_best = new_ffunc_value;
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

vector<int> IteratedLocalSearchParallel(string tai_filename, int iterations, int pert_coef,
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
	location2factory = BestLocalSearchParallel(location2factory, data, fitness_function_value);
	vector<int> location2factory_best = location2factory;

	for (int i = 0; i < iterations; ++i) {
		// perturbation
		perturbation(location2factory, pert_coef);

		// local search
		long long fitness_function_value_perturb = 0;
		vector<int>& location2factory_perturb = BestLocalSearchParallel(location2factory, data,
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

static vector<int> makeRandomStartBatch(const int number_of_ver){
    vector<int> range(number_of_ver);
    vector<int> result(number_of_ver);
    for(int i = 0; i < number_of_ver; ++i) {
        range[i] = i;
    }
    for(int i = 0; i < number_of_ver; ++i) {
            auto el = rand() % (number_of_ver - i);
            result[i]=range[el];
            swap(range[el], range[number_of_ver - i - 1]);
    }
    return result;
};

static void twoOptSwap(vector<int>& vec){
    int first = rand() % vec.size();
    int second = rand() % vec.size();
    if (second == first && first < vec.size()-1) first++;
    else first--;
    swap(vec[first], vec[second]);
}


//first-improvement
vector<int> SimpleLocalSearch(string tai_filename, int iterations,
                              long long& best_value){
    // getting the data
    QAP_data data(tai_filename);
    int size = data.size();
    // generating a random solution
    bool flag = false;
    auto sol = makeRandomStartBatch(size);
    long long val = data.fitness_function(sol);
   // vector<int> dont_look_bits(size, 0);
    for(size_t n = 0; n < iterations; ++n){
        vector<int> dont_look_bits(size, 0);
        for(size_t i = 0; i < size; ++i) {
            if(flag) {
                flag = false;
                break;
            }
            for (size_t j = 0; j < size; ++j) {
                if(i!=j && dont_look_bits[j]==0) {
                  //  auto tmp_val = data.recalculateFitness(sol, i, j,  val);
                    swap(sol[i], sol[j]);
                    long long tmp_val = data.fitness_function(sol);
                    if (tmp_val < val) {
                        val = tmp_val;
                        flag = true;
                        break;
                    }
                    else
                        swap(sol[j], sol[i]);
                }
                if(j==size-1){
                    dont_look_bits[i]=1;
                }
            }
        }
    }
    best_value = val;
    return sol;
}
