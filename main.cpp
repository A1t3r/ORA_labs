#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <queue>
#include <cmath>
#include <math.h>
#include <fstream>


struct Node {
	map<char, Node*> links;
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
	vector<bool> found;
	size_t size = 0;

	void AddPattern(const string& pattern, int point) {
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

	void Print(Node* node, string prefix) {
		cout << prefix << node << "|" << node->fake_link << endl;
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
	Bor(const vector<string>& patterns) {
		size = patterns.size();
		root = new Node();

		int number = 0;
		for (const auto& pattern : patterns) {
			AddPattern(pattern, number);
			number++;
		}

		queue<Father_Son_Value*> Fathers_and_Sons;

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

	void FindIn(string& base) {
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

	vector<bool> getFoundCopy() {
		return found;
	}

	void Show() {
		Print(root, "  ");
	}

	~Bor() {
		DelNode(root);
	}
};


vector<bool> Aho_Corasic_algorithm(string& base, vector<string>& patterns) {
	Bor bor(patterns);
	bor.FindIn(base);
	return bor.getFoundCopy();
}


void run_Aho_Corasic_algorithm_demo() {
	vector<string> pats = { "acc", "ac", "cat", "gcc", "oca", "a", "tua", "tgg" };
	string base = "ocatuaccbutnotgc";
	int number = 0;
	for (auto item : Aho_Corasic_algorithm(base, pats)) {
		string buf = "not found";
		if (item == 1) buf = "found"
		cout << pats[number] << ": " << buf << endl;

		number++:
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
	run_Aho_Corasic_algorithm_demo();
	return 0;

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
