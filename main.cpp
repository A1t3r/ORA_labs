#include "AntColony.h"
#include "SymmetricMatrix.h"
#include <chrono>

using namespace std;

double get_cost(vector<int>& result, AntColonyData& ad){
    double total = 0;
    int i = 0;
    for (i = 0; i < result.size() - 1; ++i) {
        total += ad.get_road_len(result[i], result[i + 1]);
        cout << result[i] << " ";
    }
    total += ad.get_road_len(result[i], 0);
    return total;
};

int main() {
    chrono::steady_clock::time_point pr_StartTime;
    chrono::steady_clock::time_point pr_EndTime;

    vector<string> A;
    vector<string> B;
    ifstream common_file("../data/filenames.txt");
    string tmp;
    common_file >> tmp;
    common_file >> tmp;
    while (tmp != "B") {
        A.push_back(tmp);
        common_file >> tmp;
    }
    while( common_file >> tmp)
        B.push_back(tmp);

    string file_template="../data/A/";
/*
    for(const string& name : A){
        AntColonyData ad1;
        double total_time = 0;
        ifstream tmp_file(file_template + name);
        ad1.parse_data(tmp_file, name);

        pr_StartTime = std::chrono::steady_clock::now();
        for (size_t i = 0; i < 10; ++i)
            std::vector<int> result = AntAlgorithm(
                    ad1,   // data
                    1000, // iterations
                    5,    // ants num
                    1,    // elites num
                    1.0,  // alpha
                    1.3,  // beta
                    0.03  // p
            );
        pr_EndTime = std::chrono::steady_clock::now();
        total_time = (double) std::chrono::duration_cast<std::chrono::seconds>(pr_EndTime - pr_StartTime).count() / 10.;

        std::vector<int> result = AntAlgorithm(
                ad1,   // data
                1000, // iterations
                5,    // ants num
                1,    // elites num
                1.0,  // alpha
                1.3,  // beta
                0.03  // p
        );
        ad1.save_to_file(result, get_cost(result, ad1), total_time);
    }

    file_template="../data/B/";

    for(const string& name : B){
        AntColonyData ad2;
        double total_time = 0;
        ifstream tmp_file(file_template + name);
        ad2.parse_data(tmp_file, name);

        pr_StartTime = std::chrono::steady_clock::now();
        for (size_t i = 0; i < 10; ++i)
            std::vector<int> result = AntAlgorithm(
                    ad2,   // data
                    1000, // iterations
                    5,    // ants num
                    1,    // elites num
                    1.0,  // alpha
                    1.3,  // beta
                    0.03  // p
            );
        pr_EndTime = std::chrono::steady_clock::now();
        total_time = (double) std::chrono::duration_cast<std::chrono::seconds>(pr_EndTime - pr_StartTime).count() / 10.;

        std::vector<int> result = AntAlgorithm(
                ad2,   // data
                1000, // iterations
                5,    // ants num
                1,    // elites num
                1.0,  // alpha
                1.3,  // beta
                0.03  // p
        );
        ad2.save_to_file(result, get_cost(result, ad2), total_time);
    }
*/
    SymmetricMatrix sm;
    ifstream tmp_file("../data/A/A-n32-k5.vrp");
    AntColonyData ad;
    ad.parse_data(tmp_file, "A-n32-k5.vrp");
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

