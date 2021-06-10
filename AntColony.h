
#ifndef ORA_LABS_ANTCOLONY_H
#define ORA_LABS_ANTCOLONY_H

#include "SymmetricMatrix.h"

#include <string>
#include <vector>


static double get_len(std::pair<int, int>& x, std::pair<int, int>& y){
    return sqrt(pow(std::get<0>(y) - std::get<0>(x), 2) +
                                pow(std::get<1>(y) - std::get<1>(x), 2));
}

class AntColonyData {
private:
    SymmetricMatrix SM;
    int dimension;
    int capacity;
    size_t start_point_it = 0;
    std::vector<int> demand;
public:
    AntColonyData(){};
    ~AntColonyData(){};
    void parse_data(std::ifstream& fin);

};

void AntColonyData::parse_data(std::ifstream& fin) {
    std::string tmps;
    // get common information
    while(tmps!="NODE_COORD_SECTION") {
        fin >> tmps;
        if (tmps == "DIMENSION"){
            fin >> tmps;;
            fin >> this->dimension;
        }
        if (tmps == "CAPACITY"){
            fin >> tmps;;
            fin >> this->capacity;
        }
    }
    SM.initialize(dimension);
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
                SM.item(i, j) = get_len(tmp_coords[i], tmp_coords[j]);
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
    }
}


#endif //ORA_LABS_ANTCOLONY_H
