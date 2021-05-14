#include <iostream>

#include "Object.h"
#include "BinaryChromosome.h"
#include <string>
#include <fstream>
#include <chrono>

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

double salesman_problem_genalg(vector<pair<double, double>>& map, int population_size, int it_population_size, size_t iterations){
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
    for(auto item : vecres) {
        cout << item << " ";
    }
    cout<<endl;
    return res;
}

double salesman_problem_genalg(vector<vector<double>>& map, int population_size, int it_population_size, size_t iterations){
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
    for(auto item : vecres) {
        cout << item << " ";
    }
    cout<<endl;
    return res;
}

int main() {
 //   vector<pair<int, int>> data= {{0, 0}, {0, 5}, {2, 3}, {5, 5}, {8, 3}, {6, 0}};
//    vector<vector<double>> data= {{0, 34, 36, 37, 31, 33, 35},
 //                              {34, 0, 29, 21, 22, 25, 24},
   //                            {36, 29, 0, 17, 12, 18, 17},
  //                             {37, 23, 17, 0, 32, 30, 29},
  //                             {31, 22, 12, 32, 0, 26, 24},
  //                             {33, 25, 18, 30, 26, 0, 19},
 //                              {35, 24, 17, 29, 24, 19, 0}};
 //   auto check = make_start_batch(10, 5, 0);
  //  auto check2 = get_len({0,0},{4,3});
//    auto check3 = TS_selection(check, data, 7);
//    vector<int> c1 = {8,2,3,7,1,6,0,5,4,8};
//    vector<int> c2 = {8,3,1,4,0,5,7,2,6,8};
//    partially_mapped_crossover(c1, c2);
 //   cout<<salesman_problem_genalg(data, 7, 7, 50)<<endl;
//    return 0;
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
        pr_StartTime = std::chrono::steady_clock::now();
        cout<<"answer - "<<salesman_problem_genalg(map, map.size(), map.size()/2, 70)<<endl;
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count();
        cout<<endl;
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
        pr_StartTime = std::chrono::steady_clock::now();
        cout<<"answer - "<<salesman_problem_genalg(map, map.size(), map.size()/2, 70)<<endl;
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count();
        cout<<endl;
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
        pr_StartTime = std::chrono::steady_clock::now();
        cout<<"answer - "<<salesman_problem_genalg(map, map.size(), map.size()/2, 70)<<endl;
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count();
        cout<<endl;
    }

	//knapsack_genalg(example, 10);
	return 0;
}
