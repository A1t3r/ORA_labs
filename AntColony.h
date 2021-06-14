
#ifndef ORA_LABS_ANTCOLONY_H
#define ORA_LABS_ANTCOLONY_H

#include "SymmetricMatrix.h"

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <set>


static double get_len(std::pair<int, int>& x, std::pair<int, int>& y){
    return sqrt(pow(std::get<0>(y) - std::get<0>(x), 2) +
                                pow(std::get<1>(y) - std::get<1>(x), 2));
}

class Ant{
    std::vector<int> path;
    std::set<int> unvisited;
    double total_len;
    int capacity;
    int curr_point;
public:
    void set_unvisited(int points_num) {
        for (int i = 0; i < points_num; ++i) {
            unvisited.insert(i);
        }
    }

    std::vector<int>& get_path() {
        return path;
    }

    double get_total_len() {
        return total_len;
    }

    const std::set<int>& get_unvisited() {
        return unvisited;
    }

    void add_point(int point, int capacity, double path_len) {
        path.push_back(point);
        this->capacity -= capacity;
        total_len += path_len;
        unvisited.erase(point);
        curr_point = point;
    }

    void add_start_point(int point, int capacity, double path_len) {
        path.push_back(point);
        this->capacity -= capacity;
        total_len += path_len;
        curr_point = point;
    }

    int get_curr_point() {
        return curr_point;
    }

    void set_capacity(int cap){
        capacity=cap;
    }

    int get_capacity() {
        return capacity;
    }
};

class AntColonyData {
private:
    SymmetricMatrix roads_len;
    SymmetricMatrix pheromone_index;
    int dimension;
    int capacity;
    size_t start_point_it = 0;
    std::vector<int> demand;
    float alpha, beta;
    double overall_need = 0;
    double optim = -1;
    int number_of_trucks;
public:
    AntColonyData(){};
    ~AntColonyData(){};
    void parse_data(std::ifstream& fin);
    void initilaze_pheromone();
    double get_road_len(size_t i, size_t j);
    double get_pheromone(size_t i, size_t j);
    int get_demand(size_t i);
    int get_capacity();
    void set_params(float alpha, float beta);
    int get_dimension();
    size_t get_start_point();
    int get_number_of_trucks();
    void save_to_file(std::vector<std::vector<int>>){};
    void fade_all_pheromones(double p = 0.05);
    void update_pheromone(int i, int j, double p = 0.05);
};

int AntColonyData::get_number_of_trucks(){
    return number_of_trucks;
};

int AntColonyData::get_demand(size_t i){
    return demand[i];
};

int AntColonyData::get_capacity(){
    return capacity;
};

size_t AntColonyData::get_start_point(){
    return start_point_it;
};

int AntColonyData::get_dimension(){
    return this->dimension;
};

void AntColonyData::set_params(float alpha, float beta){
    this->alpha = alpha;
    this->beta = beta;
};

double AntColonyData::get_road_len(size_t i, size_t j){
    return roads_len.item(i, j);
};

double AntColonyData::get_pheromone(size_t i, size_t j){
    return pheromone_index.item(i, j);
};

void AntColonyData::initilaze_pheromone(){
    pheromone_index.initialize(this->dimension);
    for(int i =0; i<dimension; ++i){
        for(int j =0; j<dimension-i; ++j) {
            pheromone_index.item(i,j) = (double) 1 / dimension;
        }
    }
}

void AntColonyData::fade_all_pheromones(double p) {
    double fading = 1 - p;
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension - i; ++j) {
            pheromone_index.item(i, j) *= fading;
        }
    }
}

void AntColonyData::update_pheromone(int i, int j, double p) {
    pheromone_index.item(i, j) += p * (double) 1 / dimension;
}

void AntColonyData::parse_data(std::ifstream& fin) {
    std::string tmps;
    // get common information
    while(tmps!="NODE_COORD_SECTION") {
        fin >> tmps;
        if (tmps == "trucks:"){
            fin >> tmps;
            std::string tmps2 = tmps.substr(0, tmps.size()-1);
            this->number_of_trucks = std::stoi(tmps);
        }
        else if (tmps == "value:"){
            fin >> tmps;
            std::string tmps2 = tmps.substr(0, tmps.size()-1);
            this->optim = std::stoi(tmps);
        }
       else if (tmps == "DIMENSION"){
            fin >> tmps;;
            fin >> this->dimension;
        }
       else if (tmps == "CAPACITY"){
            fin >> tmps;;
            fin >> this->capacity;
        }
    }
    roads_len.initialize(dimension);
    int value = 0;
    std::pair<int, int> first_point;
    std::pair<int, int> second_point;
    // getting len between points
    std::vector<std::pair<int, int>> tmp_coords(dimension);
    for (int i = 0; i < dimension; ++i) {
        fin >> value;
        fin >> first_point.first;
        fin >> first_point.second;
        tmp_coords[value-1]=first_point;
        }

    for (int i = 0; i < tmp_coords.size(); ++i) {
        for (int j = 0; j < tmp_coords.size(); ++j) {
             if(i <= j)
                 roads_len.item(i, j) = get_len(tmp_coords[i], tmp_coords[j]);
        }
    }

    // getting demand section
    demand.resize(dimension);
    while(tmps!="DEMAND_SECTION"){
        fin>>tmps;
    }
    int iter = 0;
    for (int i = 0; i < dimension; ++i) {
        fin >> iter;
        fin >> value;
        if(value == 0) this->start_point_it = iter-1;
        demand[iter-1]=value;
        overall_need += value;
    }
}


// ALGORITHM

int GetNextPoint(AntColonyData& data, Ant& ant, double alpha, double beta) {

    std::vector<double> P_upper_coefs;
    std::vector<int> unvisited_available;
    double P_upper_sum = 0;
    double curr_value = 0;
    int curr_point = ant.get_curr_point();
    int ant_capacity = ant.get_capacity();
    
    // choose where demand <= ant_capacity
    for (const int& unvisited : ant.get_unvisited()) {
        if (data.get_demand(unvisited) <= ant_capacity) {
            unvisited_available.push_back(unvisited);
        }
    }

    // not found any available
    if (unvisited_available.size() == 0) return -1;

    for (int& unvisited : unvisited_available) {
        curr_value = pow(data.get_pheromone(curr_point, unvisited), alpha) *
            pow((1 / data.get_road_len(curr_point, unvisited)), beta);

        P_upper_sum += curr_value;
        P_upper_coefs.push_back(curr_value);
    }

    const double accurasy = 100.;

    double mult_coef = accurasy / P_upper_sum;

    double choose_value = rand() % (int)floor(P_upper_sum * mult_coef);
    int index = -1;
    while (choose_value >= 0) {
        ++index;
        choose_value -= P_upper_coefs[index] * mult_coef; // might be out of range due to accum error
    }

    return unvisited_available[index];
}

void GlobalUpdate(AntColonyData& data, std::vector<Ant>& ants, int number_of_elites, double p) {
    sort(ants.begin(), ants.end(), [](Ant& ant_l, Ant& ant_r) -> bool {
        return ant_l.get_total_len() > ant_r.get_total_len();
        }
    );

    if (number_of_elites > ants.size()) number_of_elites = ants.size();

    // update number_of_elites best ants
    for (int i = 0; i < number_of_elites; ++i) {
        auto path = ants[i].get_path();

        for (int j = 0; j < path.size() - 1; ++j) {
            data.update_pheromone(path[j], path[j+1], p);
        }
    }
}

std::vector<int> AntAlgorithm(AntColonyData& data, int number_of_its, int number_of_ants, int number_of_elites,
    double alpha, double beta, double p) {

    srand(time(0));

    double best_found_value = -1;
    std::vector<int> best_path;

    data.initilaze_pheromone(); // init with 1/dimension probability
    int dim = data.get_dimension();

    for (size_t it = 0; it < number_of_its; ++it) {
        int start_point = data.get_start_point();

        std::vector<Ant> ants(number_of_ants);

        // fading pheromones
        data.fade_all_pheromones();

        for (Ant& ant : ants) {
            // just acceleration for the next point search
            ant.set_unvisited(dim);

            ant.set_capacity(data.get_capacity());
            ant.add_point(start_point, 0, 0);

            // choosing next point with available check
            while (ant.get_unvisited().size() != 0) {
                int next_point = GetNextPoint(data, ant, alpha, beta);

                int curr_point = ant.get_curr_point();

                if (next_point == -1) { // no availabe point to go
                    next_point = 0;
                    ant.add_start_point(
                        start_point, 0,
                        data.get_road_len(ant.get_curr_point(), start_point)
                    );
                    ant.set_capacity(data.get_capacity());
                }
                else {
                    ant.add_point(next_point, data.get_demand(next_point),
                        data.get_road_len(ant.get_curr_point(), next_point));
                }

                // local pheromone update
                data.update_pheromone(curr_point, next_point, p);
            } // one ant finished

            // get back to the start
            ant.add_start_point(
                start_point, 0,
                data.get_road_len(ant.get_curr_point(), start_point)
            );

            // choose best from best and current
            double ant_result = ant.get_total_len();
            if (ant_result < best_found_value || best_found_value == -1) {
                best_found_value = ant_result;
                best_path = ant.get_path();
            }
        } // all the ants finished

        // global pheromone update
        GlobalUpdate(data, ants, number_of_elites, p);
    }

    std::cout << best_found_value << std::endl;
    return best_path;
}

#endif //ORA_LABS_ANTCOLONY_H
