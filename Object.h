#pragma once
#include <vector>


struct Object {
	int id = 0;
	int cost = 0;
	int weight = 0;
	double quality = 0.0;
	float used_scale = 0.0;
	bool used = false;
	short use_condition = -1;
	Object(int obj_cost, int obj_weight, int obj_id) :
		cost(obj_cost), weight(obj_weight), id(obj_id) {
	};
	Object(int obj_cost, int obj_weight) :
		cost(obj_cost), weight(obj_weight) {
	};
};


std::vector<Object> example = {
        { 4 , 1 },
        { 3 , 2 },
        { 2 , 6 },
        { 4 , 2 },
};