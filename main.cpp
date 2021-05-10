#include <iostream>

#include "Object.h"
#include "BinaryChromosome.h"
#include <string>

using namespace std;

void show_population(vector<BinaryChromosome> population, string message="") {
	return;
	cout << message << endl;
	for (auto& item : population) {
		cout << item << endl;
	}
}

int knapsack_fitness_function(BinaryChromosome& bc, vector<Object> objs, int W) {
	int value = 0;
	int weight = 0;
	for (int i = 0; i < objs.size(); ++i) {
		if (bc[i]) {
			value += objs[i].cost;
			weight += objs[i].weight;
		}
	}
	if (weight > W)
		return 0;
	return value;
}

void RouletteWheelSelection(vector<BinaryChromosome>& population, 
	vector<int>& ffvalues, int ffvalues_sum, int new_population_size) {
	vector<BinaryChromosome> new_population = {};
	new_population.reserve(population.size());
	int chromosome_index = 0;

	for (int i = 0; i < new_population_size; i++) {
		int rand_int = rand() % ffvalues_sum + 1;
		chromosome_index = 0;
		while (true) {
			if (chromosome_index > ffvalues.size()) throw -1;

			rand_int -= ffvalues[chromosome_index];
			if (rand_int <= 0) {
				break;
			}
			++chromosome_index;
		}
		new_population.push_back(population[chromosome_index]);
	}

	swap(new_population, population);
	return;
}

vector<bool>& knapsack_genalg(const vector<Object>& objects, int W, int iterations=100,
	int population_size=10, int new_population_size = 6, float mutation_chance=0.01) {
	vector<bool> ret = {};

	int best_cost = -1;
	BinaryChromosome best_chromosome;

	// checks
	if (objects.empty()) return ret;
	if (new_population_size < 2) throw -2;

	// creating population
	BinaryChromosome::size = objects.size();
	vector<BinaryChromosome> population(new_population_size);
	show_population(population, "created");

	for (int i = 0; i < iterations; ++i) {
		// crossover
		for (int j = 0; j < (population_size - new_population_size + 1) / 2; ++j) {
			int first_parent = rand() % new_population_size;
			int second_parent = 0;
			do {
				second_parent = rand() % new_population_size;
			} while (second_parent == first_parent);
			population.push_back(population[first_parent]);
			population.push_back(population[second_parent]);
			population[new_population_size + j * 2].OnePointCrossover(population[new_population_size + j * 2 + 1]);
		}
		show_population(population, "crossover");

		// mutation 
		for (auto& chromosome : population) {
			chromosome.SimpleMutation(mutation_chance);
		}
		show_population(population, "mutation");

		// selection
		vector<int> costs(population_size);
		int cost_sum = 0;

		for (int j = 0; j < population_size; ++j) {
			costs[j] = knapsack_fitness_function(population[j], objects, W);
			cost_sum += costs[j];

			if (costs[j] > best_cost) {
				best_chromosome = population[j];
				best_cost = costs[j];
			}
		}

		RouletteWheelSelection(population, costs, cost_sum, new_population_size);
		show_population(population, "selected");
	}
	cout << "Chromosome: " << best_chromosome << " with cost " << best_cost << endl;
	return best_chromosome.GetGenes();
}

int main() {
	knapsack_genalg(example, 10);
	return 0;
}
