#include "AntColony.h"
#include "SymmetricMatrix.h"

using namespace std;

int main(){
    SymmetricMatrix sm;
    ifstream tmp_file("../data/A/A-n32-k5.vrp");
    AntColonyData ad;
    ad.parse_data(tmp_file);
    std::vector<int> result = AntAlgorithm(
        ad,   // data
        1000, // iterations
        5,    // ants num
        1,    // elites num
        1.0,  // alpha
        1.3,  // beta
        0.03  // p
    );

    double total = 0;
    int i = 0;
    for (i = 0; i < result.size() - 1; ++i) {
        total += ad.get_road_len(result[i], result[i + 1]);
        cout << result[i] << " ";
    }
    total += ad.get_road_len(result[i], 0);
    cout << result[i] << " ";

    cout << endl << total << endl;


    return 0;
}

