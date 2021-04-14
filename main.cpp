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
        if (fsv->father == root) {
            fsv->son->fake_link = root;
            return;
        }

        Node* father = fsv->father->fake_link;

        while (father->links.find(fsv->value) == father->links.end()) {
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
                if (current_node->links.find(item) != current_node->links.end()) {
                    current_node = current_node->links[item];
                    break;
                }
                else {
                    if (current_node == root) {
                        break;
                    }
                    current_node = current_node->fake_link;
                }
            }

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

int get_pol_hash(std::string& base, size_t len, int p){
    int hash = 0;
    for(size_t i = 0; i < len; ++i) {
        hash += (int)base[i] * pow(p, i);
    }
    return hash;
}

int recalculate_hash(std::string& base, size_t start, size_t mask_size, int prev_hash, int p){
    return (prev_hash-(int)base[start-1])/p + base[start+mask_size-1]*pow(p,mask_size-1);
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
    int p = 901;
    int mask_hash = get_pol_hash(mask, mask.size(), p)%p;
    std::vector<int> table_of_hashs{get_pol_hash(text, mask.size(), p)};
    for(size_t i = 1; i<text.size()-mask.size()+1; ++i){
        table_of_hashs.push_back(recalculate_hash(text, i, mask.size(), table_of_hashs[i-1],p));
    }
    for(size_t i = 0; i<table_of_hashs.size(); ++i){
        if(table_of_hashs[i]%p == mask_hash)
            if(check_string_part(text, mask, i))
                return i;
    }
    return -1;
}

int main() {
//	run_Aho_Corasic_algorithm_demo();
//	return 0;
    setlocale(LC_ALL, "rus");
    std::string text;
    std::string mask="ак";
   // std::cout<<mask<<std::endl;
   // std::cout<< (unsigned int)static_cast<unsigned char>(mask[0])<<std::endl;
   // return 0;

    text = "abcabaabcabca";
    mask = "abaa";
    std::cout << naive_search(text, mask) << std::endl;
    std::cout << AMBH(text, mask) << std::endl;
    std::cout << ARK(text, mask) << std::endl;
    text = "personal daata";
    mask = "daata";
    std::cout << naive_search(text, mask) << std::endl;
    std::cout << AMBH(text, mask) << std::endl;
    std::cout << ARK(text, mask) << std::endl;

    std::chrono::steady_clock::time_point pr_StartTime;
    std::chrono::steady_clock::time_point pr_EndTime;

    int (*func_arr[3])(std::string&, std::string&)={naive_search, AMBH, ARK};
    std::vector<std::string>names={"naive","AMBH", "ARK"};

    std::string filename_text = "../benchmarks/bad_t_";
    std::string filename_template = "../benchmarks/bad_w_";
    std::cout<<"Now bad template and bad text are testing..."<<std::endl;

    for(int n = 0; n < 3; ++n){
        std::cout << names[n] << " - is now checking" << std::endl;
        for(int j = 1; j < 5; ++j) {
            std::cout << "now testing sample with number " << j << "\n";
            std::ifstream file_text(filename_text + std::to_string(j) + ".txt");
            std::ifstream file_template(filename_template + std::to_string(j) + ".txt");
            file_text >> text, file_template >> mask;
            std::cout<<func_arr[n](text, mask)<<"\n";
            pr_StartTime = std::chrono::steady_clock::now();
            for (int i = 0; i < 10; ++i)
                func_arr[n](text, mask);
            pr_EndTime = std::chrono::steady_clock::now();
            std::cout << " total time is = "
                      << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10
                      << std::endl;
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
                  << std::endl;
    }
    std::string tmp_text;
    std::string tmp_template;
    filename_text = "../benchmarks/good_t_";
    filename_template = "../benchmarks/good_w_";
    std::cout<<"\n\n\nNow good template and good text are testing..."<<std::endl;
    std::ifstream file_text (filename_text + "4.txt");
    std::ifstream file_template (filename_template + "4.txt");
    for(int n = 0; n < 3; ++n){
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
       //     std::cout<<text<<std::endl;
       //     std::cout<<mask<<std::endl;
            std::cout<<func_arr[n](text, mask)<<"\n";
            pr_StartTime = std::chrono::steady_clock::now();
            for (int i = 0; i < 10; ++i)
                func_arr[n](text, mask);
            pr_EndTime = std::chrono::steady_clock::now();
            std::cout << " total time is = "
                      << std::chrono::duration_cast<std::chrono::microseconds>(pr_EndTime - pr_StartTime).count() / 10
                      << std::endl;
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
                  << std::endl;
    }
    return 0;

}
