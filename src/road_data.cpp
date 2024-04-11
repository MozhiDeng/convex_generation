#include"road_data.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

//Load data
std::vector<std::vector<float>> load_data()
{
	std::vector<std::vector<float>> data;
	std::ifstream f;
	f.open("/home/dmz/Baidu/convex_generation/data1.txt");
	std::string str;
	while (getline(f, str))			
	{
		std::istringstream ss(str);  
		std::vector<float> tmp;
		float a;
		while (ss)
		{
			std::string s;
			if (!getline(ss, s, ','))
			{
				break;
			}
			a = stof(s);
			tmp.push_back(a);
		}
		data.push_back(tmp);
	}
	return data;
}

//Print road data information
void print_vector(std::vector<std::vector<float>>& data)
{
	std::vector<float>::iterator it;
	std::vector<std::vector<float>>::iterator iter;
	std::vector<float>vec_tmp;

	for (iter = data.begin(); iter != data.end(); iter++)
	{
		vec_tmp = *iter;
		for (it = vec_tmp.begin(); it != vec_tmp.end(); it++)
		{
			std::cout << *it << " ";
		}
		std::cout << "\n";
	}
}