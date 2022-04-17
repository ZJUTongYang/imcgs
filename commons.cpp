#include "commons.h"
#include <vector>
#include "graph_problem.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <iterator>

int bestColor(const int* adj_cell, const int adj_cell_num, int* color, const int color_num, const int* cur_state)
{
	long long int bit_color = 0;
	const long long int bit_1 = 1;
	for(unsigned int i = 0; i < color_num; ++i)
	{
		bit_color |= bit_1 << color[i];
	}
	
	long long int bit_choice = 0;
	for(unsigned int i = 0; i < adj_cell_num; ++i)
	{
		if (!(bit_choice & bit_1 << cur_state[adj_cell[i]]))
		{
			bit_choice |= bit_1 << cur_state[adj_cell[i]];
		}
		else
		{
			if (bit_color & bit_1 << cur_state[adj_cell[i]])
			{
	
				return cur_state[adj_cell[i]];
			}
		}
	}
	
	if (!(bit_color & bit_choice))
	{
		return color[0];
	}
	
	// We cannot find a good colour to reduce the number of 
	for(unsigned int i = 0; i < color_num; ++i)
	{
		if (bit_choice & bit_1 << color[i])
			return color[i];
	}
}

int cal_num_of_used_colours(const int* S, int num_of_cell)
{
	// We assume that the indices of colours are not larger than 64
	long long int bit_colour = 0;
	const long long int bit_1 = 1;
	int cost = 0;
	for (unsigned int i = 0; i < num_of_cell; ++i)
	{
		if (!(bit_colour & (bit_1 << S[i])))
		{
			bit_colour |= bit_1 << S[i];
			cost += 1;
		}
	}
	return cost;
}


void CellPartition(std::vector<std::vector<int> >& result, const int* preserve, const int* all_adj_cell, const int edge_num)
{
	std::vector<int> indices;
	indices.resize(edge_num, 1);

	// Each unconnectable cells (location in the adjacent cell sequences) are marked as -1
	for (unsigned int i = 0; i < edge_num; ++i)
	{
		indices[i] = preserve[i];
	}

	// We find the first partitionable adjacent cell to start with
	int loc = -1;
	for (int i = 0; i < indices.size(); ++i)
	{
		if (indices[i] != -1)
		{
			loc = i;
			break;
		}
	}

	// If we cannot find any adjacent cell to be partitioned, we return
	int thisstart = loc;
	if (loc == -1)
	{
		return;
	}

	bool partition1 = (thisstart == 0);

	for (int i = thisstart + 1; i < indices.size(); ++i)
	{
		if (thisstart == -1 && indices[i] != -1)
		{
			thisstart = i;
			continue;
		}
		if (thisstart != -1 && indices[i] != 1)
		{
			std::vector<int> temp;
			temp = linspace(thisstart, i - 1);
			result.emplace_back(temp);
			thisstart = -1;
		}
	}

	std::vector<int> temp1 = linspace(thisstart, edge_num - 1);

	if (partition1 && thisstart != -1)
	{
		result[0].insert(result[0].begin(), temp1.begin(), temp1.end());
	}
	else if (thisstart != -1)
	{
		result.emplace_back(temp1);
	}
}

int colorIntersect(int* set1, const int set1_num, int* set2, const int set2_num, int* result)
{
	if (set1_num == 0 || set2_num == 0)
		return 0;

	long long int bit_colour = 0;
	const long long int bit_1 = 1;
	for (unsigned int i = 0; i < set1_num; ++i)
	{
		bit_colour |= (bit_1 << set1[i]);
	}

	int result_num = 0;
	for (unsigned int i = 0; i < set2_num; ++i)
	{
		if (bit_colour & (bit_1 << set2[i]))
		{
			result[result_num] = set2[i];
			result_num++;
		}
	}

	return result_num;
}

void coutAVector(const int* the, int n)
{
	std::cout << "[";
	for (unsigned int i = 0; i < n; ++i)
	{
		std::cout << the[i] << ", ";
	}
	std::cout << "]" << std::endl;
}

void coutAVector(const std::vector<int>& the)
{
	std::cout << "[";
	for (auto iter = the.begin(); iter != the.end(); ++iter)
	{
		std::cout << *iter << ", ";
	}
	std::cout << "]" << std::endl;
}

void coutProblem(GraphProblem* P)
{
	for(unsigned int i = 0; i < P->num_all_cell_; ++i)
	{
		std::cout << "cell " << i << ": edge_num " << P->cell_[i].num_edge_ << ", color_num " << P->cell_[i].num_color_ << ": [";
		for(unsigned int j = 0; j < P->cell_[i].num_color_; ++j)
		{
			std::cout << P->cell_[i].possible_color_[j] << ", ";
		}
		std::cout << "], adj_edge: [";
		for(unsigned int j = 0; j < P->cell_[i].num_edge_; ++j)
		{
			std::cout << P->cell_[i].connect_edge_[j] << ", ";
		}
		std::cout << "], isRemoved: " << P->cell_[i].isRemoved_ << std::endl;
	}
	std::cout << "changes: " << P->changes_ << std::endl;
	std::cout << "cell_process: " << P->cell_subdivision_process_ << std::endl;
	std::cout << "Edges: num_all_edge_ " << P->num_all_edge_ << std::endl;
	for(unsigned int i = 0; i < P->num_all_edge_; ++i)
	{
		std::cout << i << ": " << "[" << P->edge_[i*3] << ", " << P->edge_[i*3+1] << ", " << P->edge_[i*3+2] << "]," << std::endl;
	}
	std::cout << std::endl;
}

GraphProblem* demo_read_from_file(std::string file_name)
{
	std::ifstream infile;
	std::cout << file_name.data() << std::endl;
	infile.open(file_name.data(), std::ios::in);
	assert(infile.is_open());

	std::string s;
	std::getline(infile, s);
	std::istringstream is(s);


	int cell_num;
	int edge_num;
	is >> cell_num;
	is >> edge_num;

	GraphProblem* pResult = new GraphProblem;

	pResult->num_all_cell_ = cell_num;
	pResult->cell_ = new Cell[cell_num];
	for (unsigned int i = 0; i < cell_num; ++i)
	{
		Cell& this_cell = pResult->cell_[i];

		{
			std::vector<int> connect_edge_temp;
			std::string s;
			std::getline(infile, s);
			std::istringstream is(s);

			int edge_index;
			while (!is.eof())
			{
				is >> edge_index;
				connect_edge_temp.emplace_back(edge_index);
			}
			this_cell.assignEdge(&(connect_edge_temp[0]), connect_edge_temp.size());
		}

		{
			std::vector<int> possible_color_temp;
			std::string s;
			std::getline(infile, s);
			std::istringstream is(s);

			int color;
			while (!is.eof())
			{
				is >> color;
				possible_color_temp.emplace_back(color);
			}
			this_cell.assignColor(&(possible_color_temp[0]), possible_color_temp.size());
		}
		{
			//std::string s;
			//std::getline(infile, s);
			//std::istringstream is(s);
			//is >> this_cell.isRemoved_;
			this_cell.isRemoved_ = -1;
		}
	}

	pResult->num_all_edge_ = edge_num;
	for (unsigned int i = 0; i < edge_num; ++i)
	{
		{
			std::string s;
			std::getline(infile, s);
			std::istringstream is(s);

			is >> pResult->edge_[i*3];
			is >> pResult->edge_[i*3+1];
		}
	}

	pResult->cell_subdivision_process_ = 0;
	pResult->changes_.clear();

	std::cout << "Read problem from txt file finished" << std::endl;

	return pResult;
}

std::vector<int> linspace(int a, int b)
{
	std::vector<int> result;

	for (int i = a; i <= b; ++i)
	{
		result.emplace_back(i);
	}
	return result;
}


bool isEmptyColorIntersect(const int* set1, const int set1_num, const int* set2, const int set2_num)
{
	long long int bin_set1 = 0;
	const long long int bit_1 = 1;
	
	for(unsigned int i = 0; i < set1_num; ++i)
	{
		bin_set1 |= bit_1 << set1[i];
	}

	long long int bin_set2 = 0;
	for(unsigned int i = 0; i < set2_num; ++i)
	{
		bin_set2 |= bit_1 << set2[i];
	}

	long long int b = bin_set1 & bin_set2;
	return (b == 0);
}

void newDivisionGeneration(std::vector<std::vector<int> >& D, int connect_cell_num, std::vector<int> constraint)
{
	// Given the number of edges n, we find all possible connection of edges containing the first element
	if (connect_cell_num <= 3)
	{
		std::vector<int> temp;
		temp.resize(MAX_EDGE_NUM_PER_CELL, -1);
		D.clear();
		D.emplace_back(temp);
		return;

	}

	std::vector<int> connectable;
	connectable.resize(MAX_EDGE_NUM_PER_CELL, 1);

	if (!constraint.empty())
	{
		for (unsigned int i = 0; i < constraint.size(); ++i)
		{
			connectable[constraint[i]] = 0;
		}
	}

	int loc = -1;
	std::vector<int> temp;
	temp.resize(MAX_EDGE_NUM_PER_CELL, 1);
	temp[0] = loc;
	D.emplace_back(temp);
	D.reserve(1000);
	for (int i = 1; i < connect_cell_num; ++i)
	{
		if (connectable[i] != 0)
		{
			int numP = D.size();
			for (int j = 0; j < numP; ++j)
			{
				D[j][i] = 1;
				D.insert(D.end(), D[j]);
			}

			for (int j = numP; j < 2 * numP; ++j)
			{
				D[j][i] = -1;
			}

			if (i >= 2)
			{
				// We remove the unnecessary division
				// Here we only need to check from numP, because for the former numP elements a[i] must be 0
				D.erase(std::remove_if(D.begin() + numP, D.begin() + 2 * numP,
					[&](const std::vector<int>& a) {
					return (a[i] == -1 && a[i - 1] == 0 && a[i - 2] == -1);
				}), D.end());
			}

		}
	}

	D.erase(std::remove_if(D.begin(), D.end(),
		[&](const std::vector<int>& a) {
		return (a[a.size() - 1] == 0 && a[a.size() - 2] == -1);
	}), D.end());
	// Here the value is still 2. Say num = 5, then if the last element is 3 we should remove it (because conencting to adjacent cell 4 cannot singularly form a sub-cell). 

	D.erase(D.begin());
}

void switchCellInProblem(GraphProblem* P, int process_num)
{
	// We may move the three cells with the most number of edges to the end of the cell list
	if (process_num == 3)
		return;

	coutProblem(P);
	
	// Find the cell with non-1 possible colour and maximum number of adjacent cells
	int max_index = -1;
	int max_num = -1;
	for (unsigned int i = 0; i < P->num_all_cell_ - 1 - process_num; ++i)
	{
		if (P->cell_[i].num_color_ != 1)
		{
			if (P->cell_[i].num_edge_ > max_num)
			{
				max_index = i;
				max_num = P->cell_[i].num_edge_;
			}
		}
	}
	if (max_index == -1)
		return;

	// We switch cell [max_index] to the last one ([P->num_all_cell_ - 1 - process_num])
	
	int i1 = max_index;
	int i2 = P->num_all_cell_ - 1 - process_num;

	Cell cell1(P->cell_[i1]);
	Cell cell2(P->cell_[i2]);

	// We manipulate edge
	for (unsigned int i = 0; i < P->num_all_edge_; ++i)
	{
		if (P->edge_[i * 3] == i1)
		{
			P->edge_[i * 3] = i2;
		}
		else if (P->edge_[i * 3] == i2)
		{
			P->edge_[i * 3] = i1;
		}

		if (P->edge_[i * 3+1] == i1)
		{
			P->edge_[i * 3+1] = i2;
		}
		else if (P->edge_[i * 3+1] == i2)
		{
			P->edge_[i * 3+1] = i1;
		}
	}

	// We manipulate cell
	memcpy(&(P->cell_[i1]), &cell2, sizeof(Cell));
	memcpy(&(P->cell_[i2]), &cell1, sizeof(Cell));

	P->changes_ += "cell switch (";
	P->changes_ += std::to_string(i1);
	P->changes_ += "<->";
	P->changes_ += std::to_string(i2);
	P->changes_ += "),";
	
	//coutProblem(P);
	switchCellInProblem(P, ++process_num);
}
