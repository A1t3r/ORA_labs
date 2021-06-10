#ifndef ORA_LABS_SYMMETRICMATRIX_H
#define ORA_LABS_SYMMETRICMATRIX_H

#include <vector>
#include <ostream>
#include <string>
#include <algorithm>
#include <fstream>

//using namespace std;

class SymmetricMatrix
{
public:
    SymmetricMatrix();
    ~SymmetricMatrix();
    void initialize(int n);
    friend std::ostream& operator<< (std::ostream& out, SymmetricMatrix sm);
    double& item(int i, int j);
    int size();
    void read_data(std::ifstream& fin);

private:
    std::vector<double>& operator[] (const size_t idx);
    int values_size = 0;
    std::vector<std::vector<double>> values;
};

void SymmetricMatrix::initialize(int n) {
    values_size = n;
    values.reserve(n);
    for (int i = 0; i < n; ++i) {
        std::vector<double> empty(n - i, 0);
        values.push_back(empty);
    }
}

SymmetricMatrix::SymmetricMatrix() {
}


SymmetricMatrix::~SymmetricMatrix() {
}


std::ostream& operator<< (std::ostream& out, SymmetricMatrix sm) {
    int size = sm.values_size;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            int value = sm.item(i, j);
            if (value < 10 && value >= 0) {
                out << " " << value << " ";
            }
            else {
                out << value << " ";
            }
        }
        out << std::endl;
    }
    return out;
}

std::vector<double>& SymmetricMatrix::operator[] (const size_t idx){
    return values[idx];
};

double& SymmetricMatrix::item(int i, int j) {
    std::pair<int, int> min_max = std::minmax(i, j);
    return values[min_max.first][min_max.second - min_max.first];
}

int SymmetricMatrix::size() {
    return values_size;
}

void SymmetricMatrix::read_data(std::ifstream& fin) {
    int size = 32;
    std::string tmps;
    for (int i = 0; i < size; ++i) {
        fin>>tmps;
        if(i==20) fin>>size;
    }
    initialize(size);
    int value = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            fin >> value;
            if(i <= j)
                values[i][j - i] = value;
        }
    }
}

#endif //ORA_LABS_ANTCOLONY_H