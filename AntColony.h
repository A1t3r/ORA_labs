
#ifndef ORA_LABS_ANTCOLONY_H
#define ORA_LABS_ANTCOLONY_H

#include "SymmetricMatrix.h"

#include <string>
#include <vector>
#include <iostream>


static double get_len(std::pair<int, int>& x, std::pair<int, int>& y){
    return sqrt(pow(std::get<0>(y) - std::get<0>(x), 2) +
                                pow(std::get<1>(y) - std::get<1>(x), 2));
}

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
    double get_started_value();
    void save_to_file(std::vector<std::vector<int>>){};
};

int AntColonyData::get_demand(size_t i){
    return demand[i];
};
int AntColonyData::get_capacity(){
    return capacity;
};

double AntColonyData::get_started_value(){
    return 0;
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
            pheromone_index.item(i,j) = (double) 1/dimension;
        }
    }
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

double AntAlgorithm(AntColonyData& data, int number_of_its, float alpha, float beta){
    srand(time(0));
    std::vector<std::vector<int>> total_rout;
    double best_found_value = data.get_started_value();
    data.set_params(alpha, beta);
    data.initilaze_pheromone(); // init with 1/dimension probability
    double lower_sum = 0;
    double probability_sum = 1;
    auto dim = data.get_dimension();
    for(auto i =0; i<dim; ++i){
        for(auto j =0; j<dim; ++j){
            if(i!=j)
                lower_sum+=pow(data.get_pheromone(i, j), alpha) *
                        pow((double)(1/data.get_road_len(i, j)),beta);
        }
    }
    std::vector<double> probabilities(dim);
    std::vector<bool> visited(dim, false);
    for(size_t it =0; it<number_of_its; ++it){
        int ant_cap = data.get_capacity();
        for(size_t i =0; i<dim; ++i){
            for(size_t j =0; j<dim; ++j){
                if(!visited[j] && i!=j) {
                    auto tmppr = ((data.get_pheromone(i, j), alpha) *
                                  pow((double) (1 / data.get_road_len(i, j)), beta)) / lower_sum;
                    int len_to_calc = std::to_string(tmppr).size() - 1;
                    int range = pow(10, len_to_calc);
                    int tmp_var = tmppr * range;
                    if(rand()%range>tmp_var){
                        ant_cap-=data.get_demand(j);
                        visited[j] = true;
                    }
                }
            }
        }
    }

    return best_found_value;
}


#endif //ORA_LABS_ANTCOLONY_H
