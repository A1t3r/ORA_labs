#include <iostream>

#include "Object.h"
#include "BinaryChromosome.h"
#include <string>
#include <fstream>
#include <chrono>

using namespace std;

void show_population(vector<BinaryChromosome>& population, string message="") {
	cout << message << endl;
	for (auto& item : population) {
		cout << item << endl;
	}
}

int knapsack_fitness_function(BinaryChromosome& bc, vector<Object> objs, int W, 
    int& current_weight) {

    current_weight = 0;
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
    current_weight = weight;
	return value;
}

void RouletteWheelSelection(vector<BinaryChromosome>& population, 
	vector<int>& ffvalues, int ffvalues_sum, int new_population_size) {

	vector<BinaryChromosome> new_population = {};
	new_population.reserve(population.size());
	int chromosome_index = 0;

    if (ffvalues_sum == 0) {
        for (int i = 0; i < new_population_size; i++) {
            new_population.push_back(population[i]);
        }
        swap(new_population, population);
        return;
    }

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

struct KnapsackResult {
    BinaryChromosome chromosome;
    int cost;
    int weight;
    KnapsackResult(BinaryChromosome wchromosome, int wcost, int wweight) {
        chromosome = wchromosome;
        cost = wcost;
        weight = wweight;
    }
};

KnapsackResult knapsack_genalg(const vector<Object>& objects, int W, int iterations=100,
	int population_size=10, int new_population_size = 6, float mutation_chance=0.01) {
    BinaryChromosome bc(0);
    KnapsackResult ret(bc, 0, 0);

	int best_cost = -1;
    int weight = -2;
	BinaryChromosome best_chromosome;

	// checks
	if (objects.empty()) return ret;
	if (new_population_size < 2) throw -2;

	// creating population
	BinaryChromosome::size = objects.size();
	vector<BinaryChromosome> population(new_population_size);
	//show_population(population, "created");

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
		//show_population(population, "crossover");

		// mutation 
		for (auto& chromosome : population) {
			chromosome.SimpleMutation(mutation_chance);
		}
		//show_population(population, "mutation");

		// selection
		vector<int> costs(population_size);
		int cost_sum = 0;

		for (int j = 0; j < population_size; ++j) {
            int current_weight;
			costs[j] = knapsack_fitness_function(population[j], objects, W, current_weight);
			cost_sum += costs[j];

			if (costs[j] > best_cost) {
				best_chromosome = population[j];
				best_cost = costs[j];
                weight = current_weight;
			}
		}

		RouletteWheelSelection(population, costs, cost_sum, new_population_size);
		//show_population(population, "selected");
	}

	//cout << "Chromosome: " << best_chromosome << " with cost " << best_cost << endl;
    KnapsackResult result(best_chromosome, best_cost, weight);
    return result;
}

static vector<vector<int>> make_start_batch(int population_size, int number_of_ver, int start_point){
    srand(42);
    vector<int> range(number_of_ver-1);
    for(int i = 1; i < number_of_ver; ++i) {
        range[i-1] = i;
    }
    vector<vector<int>> total(population_size);
    for(size_t j = 0; j < population_size; ++j) {
        total[j].push_back(0);
        for(int i = 0; i < number_of_ver - 1; ++i) {
            auto el = rand() % (number_of_ver - i -1);
            total[j].push_back(range[el]);
            swap(range[el], range[number_of_ver - i - 2]);
        }
        total[j].push_back(0);
    }
    return total;
};

static double get_len(pair<double, double>& x, pair<double, double>& y){
    return sqrt(pow(get<0>(y) - get<0>(x), 2) + pow(get<1>(y) - get<1>(x), 2));
}


double TS_fitness_function_coord(vector<int>& path, vector<pair<double, double>>& map){
    double res = 0;
    for(size_t i = 0; i < path.size()-1; ++i){
        res += get_len(map[path[i]], map[path[i+1]]);
    }
    return res;
}

double TS_fitness_function_coord(vector<int>& path, vector<vector<double>>& map){
    double res = 0;
    for(size_t i = 0; i < path.size()-1; ++i){
        res += map[path[i]][path[i+1]];
    }
    return res;
}

static vector<vector<int>> TS_selection(vector<vector<int>>& batch, vector<vector<double>>& map, int new_population_size){
    vector<vector<int>> new_batch;
    vector<double> fit_funcs(batch.size());
    vector<pair<double, double>> range_list;
    if(new_population_size > batch.size())
        throw- 1;
    // calculate fitness
    for(size_t i = 0; i < batch.size(); ++i){
        fit_funcs[i] = TS_fitness_function_coord(batch[i], map);
        range_list.push_back({i, fit_funcs[i]});
        //    sum+=fit_funcs[i];
        //     local_max = (local_max > fit_funcs[i]) ? fit_funcs[i] : local_max;
    }
    // sort range list
    sort(range_list.begin(), range_list.end(), [](pair<double, double>& a, pair<double, double>& b) {
        return a.second < b.second;
    });

    // number of batches - 4, batch sizes 10%, 20%, 30%, 40%, probability of bathces - 50%, 30%, 15%, 5%

    for(size_t i = 0; i < new_population_size; ++i){
        int ch = rand() % 100;
        size_t its = 0;
        if(ch<50)
            its = rand() % (int)round((double)batch.size() / 10.);
        else if(ch<80)
            its = (rand() % ((int)round((double)batch.size() / 5.))) + (int)round((double)batch.size() / 10.);
        else if (ch<95)
            its = (rand() % ((int)round((double)batch.size() / 10.) * 3)) + ((int)round((double)batch.size() / 5.));
        else
            its = (rand() % ((int)round((double)batch.size() / 5.) * 2)) + ((int)round((double)batch.size() / 10.) * 3 );
        new_batch.push_back(batch[its]);
    }

    return new_batch;
}

static vector<vector<int>> TS_selection(vector<vector<int>>& batch, vector<pair<double, double>>& map, int new_population_size){
    vector<vector<int>> new_batch;
    vector<double> fit_funcs(batch.size());
    vector<pair<double, double>> range_list;
    if(new_population_size > batch.size())
        throw- 1;
    // calculate fitness
    for(size_t i = 0; i < batch.size(); ++i){
        fit_funcs[i] = TS_fitness_function_coord(batch[i], map);
        range_list.push_back({i, fit_funcs[i]});
    //    sum+=fit_funcs[i];
   //     local_max = (local_max > fit_funcs[i]) ? fit_funcs[i] : local_max;
    }
    // sort range list
    sort(range_list.begin(), range_list.end(), [](pair<double, double>& a, pair<double, double>& b) {
        return a.second < b.second;
    });

    // number of batches - 4, batch sizes 10%, 20%, 30%, 40%, probability of bathces - 50%, 30%, 15%, 5%

    for(size_t i = 0; i < new_population_size; ++i){
        int ch = rand() % 100;
        size_t its = 0;
        if(ch<50)
            its = rand() % (batch.size() / 10);
        else if(ch<80)
            its = (rand() % ((batch.size() / 5))) + (batch.size() / 10);
        else if (ch<95)
            its = (rand() % ((batch.size() / 10) * 3)) + ((batch.size() / 5));
        else
            its = (rand() % ((batch.size() / 5) * 2)) + ((batch.size() / 10) * 3 );
        new_batch.push_back(batch[its]);
    }

    return new_batch;
}

static void swap_mutation(vector<int>& c){
    if(c.size()==1) return;
    auto ch1 = (rand() % (c.size()-2)) +1;
    auto ch2 = (rand() % (c.size()-2)) +1;
    if(ch1 == ch2 && ch2!=c.size()-2) ch2 = ch1 + 1;
    else if(ch1 == ch2) ch2 = ch1 - 1;
    swap(c[ch1], c[ch2]);
}

void partially_mapped_crossover(vector<int>& chrm1, vector<int>& chrm2){
    vector<int> map(chrm1.size(), -1);
    vector<int> map_reverse(chrm1.size(), -1);
    size_t place1 = rand() % (chrm1.size()-2);
    size_t place2 = rand() % (chrm1.size()-2);
    place1++, place2++;
    if(place1 > place2){
        auto tmp = place1;
        place1 = place2;
        place2 = tmp;
    }
    for(size_t i = place1; i <= place2; ++i){
        map[chrm1[i]] = chrm2[i];
        map_reverse[chrm2[i]] = chrm1[i];
        auto tmp = chrm1[i];
        chrm1[i] = chrm2[i];
        chrm2[i] = tmp;
    }
    for(size_t i = 0; i < chrm1.size(); ++i){
        if(i<place1 || i>place2) {
            if (map_reverse[chrm1[i]] != -1) {
                auto tmp = map_reverse[chrm1[i]];
                while (map_reverse[tmp] != -1)
                    tmp = map_reverse[tmp];
                chrm1[i] = tmp;
            }
            if (map[chrm2[i]] != -1) {
                auto tmp = map[chrm2[i]];
                while (map[tmp] != -1)
                        tmp = map[tmp];
                chrm2[i] = tmp;
            }
        }
    }
}

double salesman_problem_genalg(vector<pair<double, double>>& map, int population_size,
                               int it_population_size, size_t iterations, bool show_res = false){
    vector<vector<int>> new_batch = make_start_batch(population_size, map.size(), 0);
    //vector<vector<int>> new_batch;
    for(size_t it = 0; it<iterations; ++it){
        new_batch=TS_selection(new_batch, map, it_population_size);
        for(size_t i = 0; i<population_size/2; ++i){
            auto cros_ch = rand() % 100;
            if(cros_ch < 85) {
                auto ch1 = rand() % new_batch.size();
                auto ch2 = rand() % new_batch.size();
                if (ch1 == ch2 && ch2 != new_batch.size() - 1) ch2 = ch1 + 1;
                else if (ch1 == ch2) ch2 = ch1 - 1;
                partially_mapped_crossover(new_batch[ch1], new_batch[ch2]);
            }
        }
        for(size_t i = 0; i<new_batch.size(); ++i)
            if(3 >= (rand() % 100)){
                swap_mutation(new_batch[i]);
            }
    }
    double res =TS_fitness_function_coord(new_batch[0], map);
    vector<int> vecres = new_batch[0];
    for(auto& sample : new_batch){
        double tmp = TS_fitness_function_coord(sample, map);
        if(tmp < res){
            res = tmp;
            vecres = sample;
        }
    }
    if(show_res) {
        for (auto item : vecres) {
            cout << item << " ";
        }
    }
    cout<<endl;
    return res;
}

double salesman_problem_genalg(vector<vector<double>>& map, int population_size, int it_population_size,
                               size_t iterations, bool show_res = false){
    vector<vector<int>> new_batch = make_start_batch(population_size, map.size(), 0);
    for(size_t it = 0; it<iterations; ++it){
        new_batch=TS_selection(new_batch, map, it_population_size);
        for(size_t i = 0; i<population_size/2; ++i){
            auto cros_ch = rand() % 100;
            if(cros_ch < 85) {
                auto ch1 = rand() % new_batch.size();
                auto ch2 = rand() % new_batch.size();
                if (ch1 == ch2 && ch2 != new_batch.size() - 1) ch2 = ch1 + 1;
                else if (ch1 == ch2) ch2 = ch1 - 1;
                partially_mapped_crossover(new_batch[ch1], new_batch[ch2]);
            }
        }
        for(size_t i = 0; i<new_batch.size(); ++i)
            if(3 >= (rand() % 100)){
                swap_mutation(new_batch[i]);
            }
    }
    double res =TS_fitness_function_coord(new_batch[0], map);
    vector<int> vecres = new_batch[0];
    for(auto& sample : new_batch){
        double tmp = TS_fitness_function_coord(sample, map);
        if(tmp < res){
            res = tmp;
            vecres = sample;
        }
    }
    if(show_res) {
        for (auto item : vecres) {
            cout << item << " ";
        }
    }
    cout<<endl;
    return res;
}

int knapack_main() {
    string file_template = "../knapsack_data/";
    std::chrono::steady_clock::time_point pr_StartTime;
    std::chrono::steady_clock::time_point pr_EndTime;
    int capacity = 0;
    int p = 0;
    int w = 0;
    vector<Object> objects;

    srand(0);

    // alg opts
    int iterations = 30;
    int population_size = 50; 
    int new_population_size = 20;
    float mutation_chance = 0.05;

    for (int i = 1; i < 8; ++i) {
        cout << "Now testing sample number " << i << endl;
        int id = 0;
        ifstream file_capacity(file_template + "p0" + std::to_string(i) + "_c.txt");
        file_capacity >> capacity;
        ifstream file_p(file_template + "p0" + std::to_string(i) + "_p.txt");
        ifstream file_w(file_template + "p0" + std::to_string(i) + "_w.txt");
        while (!file_p.eof()) {
            file_p >> p, file_w >> w;
            objects.emplace_back(p, w, id);
            id++;
        }

        cout << "Answer: \n";
        KnapsackResult knapsack_result = knapsack_genalg(objects, capacity, 
            iterations, population_size, new_population_size, mutation_chance);
        
        cout << knapsack_result.chromosome << ", total cost = " << knapsack_result.cost 
            << " total weight " << knapsack_result.weight << endl;

        pr_StartTime = std::chrono::steady_clock::now();
        for (size_t i = 0; i < 10; ++i)
            knapsack_genalg(objects, capacity,
                iterations, population_size, new_population_size, mutation_chance);
        pr_EndTime = std::chrono::steady_clock::now();

        std::cout << " total time is = "
            << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10
            << std::endl;

        objects.clear();
    }
    return 0;
}

int main() {
    knapack_main();

    std::chrono::steady_clock::time_point pr_StartTime;
    std::chrono::steady_clock::time_point pr_EndTime;
    string file_template = "../data/";
    vector<string> file_names = {"a280.tsp", "att48.tsp", "bays29 coord.tsp", "ch150.tsp", "fl417.tsp"};
    for (auto& file_name : file_names){
        cout<<"Now testing sample "<<file_name<<endl;
        vector<pair<double, double>> map;
        ifstream testfile(file_template+file_name);
        string tmp_text="a";
        pair<int, int> tmp_coord;
        while(!testfile.eof()) {
            getline(testfile, tmp_text);
            bool first_out = false;
            string num_s;
            for(char ch : tmp_text){
                if(ch==' '){
                    if(first_out){
                        tmp_coord.first = stod(num_s);
                    }
                    else first_out = true;
                    num_s.clear();
                }
                num_s+=ch;
            }
         //   text += tmp_text;
            tmp_coord.second = stod(num_s);
            map.push_back(tmp_coord);
        }
        cout<<"answer - "<<salesman_problem_genalg(map, map.size(), map.size()/2, 70, true)<<endl;
        pr_StartTime = std::chrono::steady_clock::now();
        for(size_t n = 0; n < 10; ++n) {
            salesman_problem_genalg(map, map.size(), map.size() / 2, 70);
        }
            pr_EndTime = std::chrono::steady_clock::now();
            std::cout << " total time is = "
                      << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10.;
            cout << endl;
        }


    vector<string> file_names2 = {"bays29.tsp"};
    for (auto& file_name : file_names2){
        cout<<"Now testing sample "<<file_name<<endl;
        vector<vector<double>> map;
        ifstream testfile(file_template+file_name);
        string tmp_text="a";
        vector<double> tmp_row;
        while(!testfile.eof()) {
            getline(testfile, tmp_text);
            string num_s;
            for(char ch : tmp_text){
                if(ch==' '){
                    tmp_row.push_back(stod(num_s));
                    num_s.clear();
                }
                num_s+=ch;
            }
            tmp_row.push_back(stod(num_s));
            map.push_back(tmp_row);
            tmp_row.clear();
        }
        cout<<"answer - "<<salesman_problem_genalg(map, map.size(), map.size()/2, 70, true)<<endl;
        pr_StartTime = std::chrono::steady_clock::now();
        for(size_t n = 0; n < 10; ++n) {
            salesman_problem_genalg(map, map.size(), map.size() / 2, 70);
        }
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10.;
        cout << endl;
    }

    vector<string> file_names3 = {"gr17.tsp"};
    for (auto& file_name : file_names3){
        cout<<"Now testing sample "<<file_name<<endl;
        vector<vector<double>> map;
        ifstream testfile(file_template+file_name);
        string tmp_text="a";
     //   vector<double> tmp_row(17, -1);
        while(!testfile.eof()) {
            int i = 0;
            vector<double> tmp_row(17, -1);
            getline(testfile, tmp_text);
            string num_s;
            for(char ch : tmp_text){
                if(ch==' '){
               //     tmp_row.push_back(stod(num_s));
                    tmp_row[i]=stod(num_s);
                    num_s.clear();
                    i++;
                }
                num_s+=ch;
              //  i++;
            }
       //     tmp_row.push_back(stod(num_s));
            tmp_row[i]=stod(num_s);
            map.push_back(tmp_row);
            tmp_row.clear();
        }
        for(size_t i = 0; i<17; ++i){
            for(size_t j = 0; j<17; ++j){
                map[i][j] = map[j][i];
            }
        }
        cout<<"answer - "<<salesman_problem_genalg(map, map.size(), map.size()/2, 70, true)<<endl;
        pr_StartTime = std::chrono::steady_clock::now();
        for(size_t n = 0; n < 10; ++n) {
            salesman_problem_genalg(map, map.size(), map.size() / 2, 70);
        }
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10.;
        cout << endl;
    }

	return 0;
}
