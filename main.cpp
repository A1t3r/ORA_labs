#include "LocalSearch.h"

#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;

int main() {
    srand(time(0));

    chrono::steady_clock::time_point pr_StartTime;
    chrono::steady_clock::time_point pr_EndTime;
    ofstream fout;
	vector<int> location2factory;
	long long best_value = -1;
	string to_files = "../data/";
   // vector<string> file_names{"tai20a","tai40a","tai60a","tai80a","tai100a"};
    vector<string> file_names{"test","test2","test3","tai20a","tai40a","tai60a","tai80a","tai100a"};
//	location2factory = IteratedLocalSearch("../data/test", 1000, 2, best_value);
//    location2factory = SimpleLocalSearch("../data/test", 100, best_value);
    for (auto& file_name : file_names) {
        cout<<"Now testing file "+file_name<<endl;

        cout<<"Iterated local search: "<<endl;
        location2factory = IteratedLocalSearch(to_files+file_name, 50, 2, best_value);
        cout << "Value " << best_value << endl << endl;
        fout.open("../out/iterated/"+file_name+".sol.txt");
        for (int i = 0; i < location2factory.size(); ++i) {
            cout << location2factory[i] << " " << i << endl;
            fout << location2factory[i]<<" ";
        }
        fout.close();
        pr_StartTime = std::chrono::steady_clock::now();
        for(size_t i =0; i < 10; ++i)
            IteratedLocalSearch(to_files+file_name, 50, 2, best_value);
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10.;
        cout << endl;
        /*
        cout<<"First-improvement local search: "<<endl;
        location2factory = SimpleLocalSearch(to_files+file_name, 200, best_value);
        cout << "Value " << best_value << endl << endl;
        fout.open("../out/first-improvement(DO NOT RUIN MY SOLUTIONS)/"+file_name+".sol.txt");
        for (int i = 0; i < location2factory.size(); ++i) {
            cout << location2factory[i] << " " << i << endl;
            fout << location2factory[i]<<" ";
        }
        fout.close();
        pr_StartTime = std::chrono::steady_clock::now();
        for(size_t i =0; i < 10; ++i)
            SimpleLocalSearch(to_files+file_name, 200, best_value);
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10.;
        cout << endl;
        */
        cout << endl;
    }

	return 0;
}
