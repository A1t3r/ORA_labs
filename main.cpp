#include <iostream>
#include <string>
#include <math.h>
#include <fstream>

static bool check_string_part(std::string& text, std::string& mask, size_t pos){
    for(size_t j = 0; j < mask.size(); ++j){
        if(text[j+pos]!=mask[j]) return false;
    }
    return true;
}

static int check_string_part_reverse(std::string& text, std::string& mask, size_t pos){
    for(int j = mask.size()-1; j >= 0; --j){
        if(text[j+pos]!=mask[j]) return j;
    }
    return -1;
}

static int get_pol_hash(std::string base, int p){
    int hash = 0;
    for(size_t i = 0; i<base.size(); ++i)
        hash += (int)base[i]*pow(p,i);
    return hash%p;
}

int naive_search(std::string& text, std::string& mask){
    for(size_t i = 0; i < text.size(); ++i){
        if (check_string_part(text, mask, i)) return i;
    }
    return -1;
}

int AMBH(std::string& text, std::string& mask){
    size_t mask_len = mask.size();
    int ASCII_table[128];
    for(size_t i = 0; i<128; ++i)
        ASCII_table[i] = mask_len;

    for(int i = mask_len-2; i>=0; --i){
        if(ASCII_table[(size_t)mask[i]]==mask_len)
            ASCII_table[(size_t)mask[i]]=mask_len-i-1;
    }

    for(size_t i = 0; i < text.size(); ++i){
        int n = check_string_part_reverse(text, mask, i);
        if(n!=-1){
            if(n+1!=mask_len)
                i+=ASCII_table[(int)mask[mask_len-1]]-1;
                //  std::cout<<text[n+i]<<std::endl;
            else
                i+=ASCII_table[(int)text[n+i]]-1;
        }
        else return i;
    }

    return -1;
}

int ARK(std::string& text, std::string& mask) {
    int p = 97;
    int mask_hash = get_pol_hash(mask, p);
//    int table_of_hashs[text.size()-mask.size()+1]; STILL IN PROGRESS
    for(size_t i = 0; i<text.size(); ++i){

        // TO DO    cutted string to pol hash
    }
    return 0;
}

int main() {
    std::string text;
    std::string mask;
    /*
    text = "abcabaabcabca";
    mask = "abaa";
    std::cout << naive_search(text, mask) << std::endl;
    std::cout << AMBH(text, mask) << std::endl;

    text = "personal daata";
    mask = "daata";
    std::cout << naive_search(text, mask) << std::endl;
    std::cout << AMBH(text, mask) << std::endl;
*/
    std::string filename_text = "../benchmarks/bad_t_";
    std::string filename_template = "../benchmarks/bad_w_";
    std::cout<<"Now bad template and bad text are testing..."<<std::endl;
    for(int n = 1; n < 5; ++n){
        std::cout<<"now testing sample with number "<<n<<"\n";
        std::ifstream file_text (filename_text+std::to_string(n) + ".txt");
        std::ifstream file_template (filename_template+std::to_string(n) + ".txt");
        file_text>>text, file_template>>mask;
        std::cout<<"results: "<<"\n";
        std::cout<<"naive "<< naive_search(text, mask) << std::endl;
        std::cout<<"AMBH "<< AMBH(text, mask) << std::endl;
    }
/*
    filename_text = "../benchmarks/good_t_";
    filename_template = "../benchmarks/good_w_";
    std::cout<<"\nNow good template and good text are testing..."<<std::endl;
    for(int n = 1; n < 5; ++n){
        std::cout<<"now testing sample with number "<<n<<"\n";
        std::ifstream file_text (filename_text+std::to_string(n) + ".txt");
        std::ifstream file_template (filename_template+std::to_string(n) + ".txt");
        file_text>>text, file_template>>mask;
        std::cout<<"results: "<<"\n";
        std::cout<<"naive "<< naive_search(text, mask) << std::endl;
        std::cout<<"AMBH "<< AMBH(text, mask) << std::endl;
    }
*/
    return 0;
}
