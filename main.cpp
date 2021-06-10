#include "AntColony.h"
#include "SymmetricMatrix.h"

using namespace std;

int main(){
    SymmetricMatrix sm;
    ifstream tmp_file("../data/A/A-n32-k5.vrp");
    AntColonyData ad;
    ad.parse_data(tmp_file);
    AntAlgorithm(ad, 5, 1, 1);
    return 0;
}

