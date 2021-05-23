#pragma once

#include <vector>
#include <ostream>
#include <string>
#include <algorithm>
#include <fstream>

using namespace std;

class SymmetricMatrix
{
public:
	SymmetricMatrix();
	~SymmetricMatrix();
	void initialize(int n);
	friend ostream& operator<< (ostream& out, SymmetricMatrix sm);
	int& item(int i, int j);
	int size();
	void read_from_taiNa(ifstream& fin, int size);

private:
    vector<int>& operator[] (const size_t idx);
	int values_size = 0;
	vector<vector<int>> values;
};

void SymmetricMatrix::initialize(int n) {
	values_size = n;
	values.reserve(n);
	//vector<int> empty(n);
	for (int i = 0; i < n; ++i) {
        vector<int> empty(n - i, 0);
		values.push_back(empty);
       // empty.reserve(n - i - 1);
	}
}

SymmetricMatrix::SymmetricMatrix() {
}


SymmetricMatrix::~SymmetricMatrix() {
}


ostream& operator<< (ostream& out, SymmetricMatrix sm) {
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
		out << endl;
	}
	return out;
}

vector<int>& SymmetricMatrix::operator[] (const size_t idx){
    return values[idx];
};

int& SymmetricMatrix::item(int i, int j) {
	pair<int, int> min_max = minmax(i, j);
	return values[min_max.first][min_max.second - min_max.first];
}

int SymmetricMatrix::size() {
	return values_size;
}

void SymmetricMatrix::read_from_taiNa(ifstream& fin, int size) {
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