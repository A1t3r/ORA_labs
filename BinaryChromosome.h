#pragma once

#include <vector>
#include <random>

using namespace std;

class BinaryChromosome {
public:
	static int size;

	BinaryChromosome();
	BinaryChromosome(const int length);
	~BinaryChromosome();

	bool&& operator[] (const int index);
	friend ostream& operator<< (std::ostream& out, const BinaryChromosome& bc);
	BinaryChromosome& operator= (const BinaryChromosome& drob);

	void OnePointCrossover(BinaryChromosome& bc);
	bool SimpleMutation(float chance);
	vector<bool>& GetGenes();

private:
	vector<bool> genes;
	const int chanse_accurasy = 1000;
};

int BinaryChromosome::size;

BinaryChromosome::BinaryChromosome() {
	genes.reserve(size);
	for (int i = 0; i < size; i++) {
		genes.push_back(rand() % 2);
	}
}

BinaryChromosome::BinaryChromosome(const int length) {
	genes.reserve(length);
	for (int i = 0; i < length; i++) {
		genes.push_back(rand() % 2);
	}
}

BinaryChromosome::~BinaryChromosome() {
}

bool&& BinaryChromosome::operator[] (const int index) {
	return genes[index];
}

void BinaryChromosome::OnePointCrossover(BinaryChromosome& bc) {
	int swap_point = rand() % genes.size();
	for (int i = swap_point; i < genes.size(); i++) {
		swap(genes[i], bc.genes[i]);
	}
	return;
}

bool BinaryChromosome::SimpleMutation(float chance) {
	bool produced_mutation = false;
	for (int i = 0; i < genes.size(); i++) {
		if ((float)(rand() % chanse_accurasy) < chance * chanse_accurasy) {
			genes[i] = (genes[i] + 1) % 2;
			produced_mutation = true;
		}
	}
	return produced_mutation;
}

ostream& operator<< (std::ostream& out, const BinaryChromosome& bc) {
	for (const bool gene : bc.genes) {
		out << gene;
	}
	return out;
}

BinaryChromosome& BinaryChromosome::operator= (const BinaryChromosome& bc) {
	this->genes = bc.genes;
	return *this;
}

vector<bool>& BinaryChromosome::GetGenes() {
	return genes;
}
