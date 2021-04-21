#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <queue>
#include <cmath>
#include <math.h>
#include <fstream>
#include <vector>
#include <chrono>
#include <locale.h>
#include <Windows.h>

size_t global_counter;

struct Node {
    std::map<char, Node*> links;
    Node* fake_link = 0;
    bool isFinish = false;
    int point = -1;
};


struct Father_Son_Value {
    Node* father;
    Node* son;
    char value;
};


class Bor {
private:
    Node* root;
    std::vector<bool> found;
    size_t size = 0;

    void AddPattern(const std::string& pattern, int point) {
        Node* current_node = root;
        Node* new_node = 0;

        for (const char& sym : pattern) {
            global_counter++;
            if (current_node->links.find(sym) == current_node->links.end()) {
                new_node = new Node();
                current_node->links[sym] = new_node;
                current_node = new_node;
            }
            else {
                current_node = current_node->links[sym];
            }
        }
        current_node->isFinish = true;
        current_node->point = point;
    }

    void DelNode(Node* node) {
        for (auto& child_node : node->links) {
            DelNode(child_node.second);
        }
        delete node;
    }

    void Print(Node* node, std::string prefix) {
        std::cout << prefix << node << "|" << node->fake_link << std::endl;
        for (auto& child_node : node->links) {
            Print(child_node.second, prefix + "  ");
        }
    }

    void SetFakeLink(Father_Son_Value* fsv) {
        global_counter++;
        if (fsv->father == root) {
            fsv->son->fake_link = root;
            return;
        }

        Node* father = fsv->father->fake_link;

        while (father->links.find(fsv->value) == father->links.end()) {
            global_counter++;
            if (father == root) {
                fsv->son->fake_link = root;
                return;
            }
            father = father->fake_link;
        }

        fsv->son->fake_link = father->links[fsv->value];
        return;
    }

public:
    Bor(const std::vector<std::string>& patterns) {
        size = patterns.size();
        root = new Node();

        int number = 0;
        for (const auto& pattern : patterns) {
            AddPattern(pattern, number);
            number++;
        }

        std::queue<Father_Son_Value*> Fathers_and_Sons;

        root->fake_link = root;
        for (auto& pair : root->links) {
            Father_Son_Value* fsv = new Father_Son_Value();
            fsv->father = root;
            fsv->son = pair.second;
            fsv->value = pair.first;
            Fathers_and_Sons.push(fsv);
        }

        while (!Fathers_and_Sons.empty()) {
            Father_Son_Value* current_fsv = Fathers_and_Sons.front();
            SetFakeLink(current_fsv);

            for (auto& pair : current_fsv->son->links) {
                Father_Son_Value* fsv = new Father_Son_Value();
                fsv->father = current_fsv->son;
                fsv->son = pair.second;
                fsv->value = pair.first;
                Fathers_and_Sons.push(fsv);
            }

            Fathers_and_Sons.pop();
            delete current_fsv;
        }

    }

    void FindIn(std::string& base) {
        for (int i = 0; i < size; i++) {
            found.push_back(false);
        }

        Node* current_node = root;

        for (auto iter = base.begin(); iter < base.end(); iter++) {
            char item = *iter;
            while (true) {
                global_counter++;
                if (current_node->links.find(item) != current_node->links.end()) {
                    current_node = current_node->links[item];
                    break;
                }
                else {
                    global_counter++;
                    if (current_node == root) {
                        break;
                    }
                    current_node = current_node->fake_link;
                }
            }
            global_counter++;
            if (current_node->isFinish) {
                found[current_node->point] = true;
            }
            Node* subnode = current_node;
            while (subnode->fake_link->isFinish) {
                found[current_node->fake_link->point] = true;
                subnode = current_node->fake_link;
            }
        }

    }

    std::vector<bool> getFoundCopy() {
        return found;
    }

    void Show() {
        Print(root, "  ");
    }

    ~Bor() {
        DelNode(root);
    }
};


std::vector<bool> Aho_Corasic_algorithm(std::string& base, std::vector<std::string>& patterns) {
    global_counter = 0;
    Bor bor(patterns);
    bor.FindIn(base);
    return bor.getFoundCopy();
}


void run_Aho_Corasic_algorithm_demo() {
    std::vector<std::string> pats = { "acc", "ac", "cat", "gcc", "oca", "a", "tua", "tgg" };
    std::string base = "ocatuaccbutnotgc";
    int number = 0;
    for (auto item : Aho_Corasic_algorithm(base, pats)) {
        std::string buf = "not found";
        if (item == 1) buf = "found";
        std::cout << pats[number] << ": " << buf << std::endl;

        number++;
    }
    return;
}

static bool check_string_part(std::string& text, std::string& mask, size_t pos){
    for(size_t j = 0; j < mask.size(); ++j){
        global_counter++;
        if(text[j+pos]!=mask[j]) return false;
    }
    return true;
}

static int check_string_part_with_pos(std::string& text, std::string& mask, size_t pos, size_t start){
    for(size_t j = start; j < mask.size(); ++j){
        global_counter++;
        if(text[j+pos]!=mask[j]) return j;
    }
    return -1;
}

static int check_string_part_reverse(std::string& text, std::string& mask, size_t pos){
    for(int j = mask.size()-1; j >= 0; --j){
        global_counter++;
        if(text[j+pos]!=mask[j]) return j;
    }
    return -1;
}

unsigned long long get_pol_hash(std::string& base, size_t len, int p){
    unsigned long long hash = 0;
    for(size_t i = 0; i < len; ++i) {
        hash += (size_t)base[i] * pow(p, i);
    }
    return hash;
}

unsigned long long get_hash(std::string& base, size_t len, int p){
    unsigned long long hash = 0;
    for(size_t i = 0; i < len; ++i) {
        hash += (size_t)base[i];
    }
    return hash;
}

unsigned long long recalculate_hash(std::string& base, size_t start, size_t mask_size, long long prev_hash, int p){
    return (prev_hash - (size_t)base[start-1])/p + (size_t)base[start+mask_size-1]*pow(p,mask_size-1);
}

unsigned long long recalculate_hash2(std::string& base, size_t start, size_t mask_size, long long prev_hash, int p) {
    return (prev_hash - (size_t)base[start - 1]) + (size_t)base[start + mask_size - 1];
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
                i+=ASCII_table[(size_t)mask[mask_len-1]]-1;
                //  std::cout<<text[n+i]<<std::endl;
            else
                i+=ASCII_table[(size_t)text[n+i]]-1;
        }
        else return i;
    }

    return -1;
}

int AMBHrus(std::string& text, std::string& mask) {
    size_t mask_len = mask.size();
    int ASCII_table[300];
    for (size_t i = 0; i < 300; ++i)
        ASCII_table[i] = mask_len;

    for (int i = mask_len - 2; i >= 0; --i) {
        if (ASCII_table[(size_t)mask[i] % 1000] == mask_len)
           // std::cout<<(size_t)mask[i] % 1000<<"\n";
            ASCII_table[(size_t)mask[i] % 1000] = mask_len - i - 1;
    }

    for (size_t i = 0; i < text.size(); ++i) {
        int n = check_string_part_reverse(text, mask, i);
        if (n != -1) {
            if (n + 1 != mask_len)
                i += ASCII_table[(size_t)mask[mask_len - 1] % 1000] - 1;
            //  std::cout<<text[n+i]<<std::endl;
            else
                i += ASCII_table[(size_t)text[n + i] % 1000] - 1;
        }
        else return i;
    }

    return -1;
}

int ARK(std::string& text, std::string& mask) {
    int p = 57;
    auto mask_hash = get_hash(mask, mask.size(), p)%p;
    std::vector<unsigned long long> table_of_hashs(text.size()-mask.size()+1);
    table_of_hashs[0]=get_hash(text, mask.size(), p);
    for(size_t i = 1; i<text.size()-mask.size()+1; ++i){
        table_of_hashs[i]=(recalculate_hash2(text, i, mask.size(), table_of_hashs[i-1], p));
    }
    for(size_t i = 0; i<table_of_hashs.size(); ++i){
        if(table_of_hashs[i]%p == mask_hash)
            if(check_string_part(text, mask, i))
                return i;
    }
    return -1;
}

std::vector<int> get_prefix_table(std::string& mask){
    std::vector<int> res(mask.size());
    res[0]=0;
    int k = 0;
    for(size_t i = 1; i<mask.size(); ++i){

        k = res[i-1];
        while(k > 0 && mask[k]!=mask[i]) {
     //       global_counter++;
            k = res[k - 1];
        }
        global_counter++;
        if(mask[k]==mask[i]) {
            res[i] = k + 1;
        }
        else res[i]=0;
    }
    return res;
}

int KMP(std::string& text, std::string& mask) {
    std::vector<int> prefix_table{get_prefix_table(mask)};
    int tmp = 0;
    int j = 0;
    for(size_t i = 0; i<text.size(); ++i) {
        global_counter++;
        if(mask[j]==text[i]){
            j++;
            if(j==mask.size()) return i - j + 1;
        }
        else if(j!=0){
            j = prefix_table[j - 1];
            i--;
        }
    }
    return -1;
}

int main() {
    setlocale(LC_ALL, "Russian");
    SetConsoleOutputCP(1251);
    SetConsoleCP(1251);

    std::string text;
    std::string mask;

    std::chrono::steady_clock::time_point pr_StartTime;
    std::chrono::steady_clock::time_point pr_EndTime;

    int (*func_arr[4])(std::string&, std::string&)={naive_search, AMBH, ARK, KMP};
    std::vector<std::string>names={"naive","AMBH", "ARK", "KMP"};

    std::string filename_text = "../benchmarks/bad_t_";
    std::string filename_template = "../benchmarks/bad_w_";
    std::cout<<"Now bad template and bad text are testing..."<<std::endl;

    for(int n = 0; n < 4; ++n){
        std::cout << names[n] << " - is now checking" << std::endl;
        for(int j = 1; j < 5; ++j) {
            std::cout << "now testing sample with number " << j << "\n";
            std::ifstream file_text(filename_text + std::to_string(j) + ".txt");
            std::ifstream file_template(filename_template + std::to_string(j) + ".txt");
            file_text >> text, file_template >> mask;
            std::cout<<func_arr[n](text, mask)<<"\n";
            pr_StartTime = std::chrono::steady_clock::now();
            global_counter = 0;
            for (int i = 0; i < 10; ++i)
                func_arr[n](text, mask);
            pr_EndTime = std::chrono::steady_clock::now();
            std::cout << " total time is = "
                      << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10
                      << " number of comparisons = " << global_counter / 10
                      << std::endl;
            global_counter = 0;
        }
    }
    std::cout  << "karasik - is now checking" << std::endl;
    for(int j = 1; j < 5; ++j) {
        std::cout << "now testing sample with number " << j << "\n";
        std::ifstream file_text(filename_text + std::to_string(j) + ".txt");
        std::ifstream file_template(filename_template + std::to_string(j) + ".txt");
        file_text >> text, file_template >> mask;
        std::vector<std::string> tmp_vec{mask};
        std::cout<<Aho_Corasic_algorithm(text, tmp_vec)[0]<<"\n";
        pr_StartTime = std::chrono::steady_clock::now();
        for (int i = 0; i < 10; ++i)
            Aho_Corasic_algorithm(text, tmp_vec);
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10
                  << " number of comparisons = " << global_counter / 10
                  << std::endl;
    }
    int (*func_arr2[4])(std::string&, std::string&) = { naive_search, AMBHrus, ARK, KMP };
    std::string tmp_text;
    std::string tmp_template;
    filename_text = "../benchmarks/goodr_t_";
    filename_template = "../benchmarks/goodr_w_";
    std::cout<<"\n\n\nNow good template and good text are testing..."<<std::endl;
    for(int n = 0; n < 4; ++n){
        std::cout << names[n] << " - is now checking" << std::endl;
        for(int j = 1; j < 5; ++j) {
            text="";
            mask="";
            std::cout << "now testing sample with number " << j << "\n";
            std::ifstream file_text(filename_text + std::to_string(j) + ".txt");
            std::ifstream file_template(filename_template + std::to_string(j) + ".txt");
            while(!file_text.eof()) {
                getline(file_text, tmp_text);
                text += tmp_text;
            }
            while(!file_template.eof()) {
                getline(file_template, tmp_template);
                mask += tmp_template;
            }
               //  std::cout<<text<<std::endl;
               //  std::cout<<mask<<std::endl;
            std::cout<<func_arr2[n](text, mask)<<"\n";
            pr_StartTime = std::chrono::steady_clock::now();
            global_counter = 0;
            for (int i = 0; i < 10; ++i)
                func_arr2[n](text, mask);
            pr_EndTime = std::chrono::steady_clock::now();
            std::cout << " total time is = "
                      << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10
                      << " number of comparisons = " << global_counter / 10
                      << std::endl;
            global_counter = 0;
        }
    }
    std::cout  << "karasik - is now checking" << std::endl;
    for(int j = 1; j < 5; ++j) {
        text="";
        mask="";
        std::cout << "now testing sample with number " << j << "\n";
        std::ifstream file_text(filename_text + std::to_string(j) + ".txt");
        std::ifstream file_template(filename_template + std::to_string(j) + ".txt");
        while(!file_text.eof()) {
            getline(file_text, tmp_text);
            text += tmp_text;
        }
        while(!file_template.eof()) {
            getline(file_template, tmp_template);
            mask += tmp_template;
        }
        std::vector<std::string> tmp_vec{mask};
        std::cout<<Aho_Corasic_algorithm(text, tmp_vec)[0]<<"\n";
        pr_StartTime = std::chrono::steady_clock::now();
        for (int i = 0; i < 10; ++i)
            Aho_Corasic_algorithm(text, tmp_vec);
        pr_EndTime = std::chrono::steady_clock::now();
        std::cout << " total time is = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10
                  << " number of comparisons = " << global_counter / 10
                  << std::endl;
    }
    return 0;

}