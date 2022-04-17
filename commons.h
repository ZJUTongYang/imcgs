#pragma once
#include "graph_problem.h"

int bestColor(const int* adj_cell, const int adj_cell_num, int* color, const int color_num, const int* cur_state);
int cal_num_of_used_colours(const int* S, int num_of_cell);

void CellPartition(std::vector<std::vector<int> >& result, const int* preserve, const int* all_adj_cell, const int edge_num);

int colorIntersect(int* set1, const int set1_num, int* set2, const int set2_num, int* result);

void coutAVector(const std::vector<int>& the);
void coutAVector(const int* the, int n);

void coutProblem(GraphProblem* P);
GraphProblem* demo_read_from_file(std::string file_name);

bool isEmptyColorIntersect(const int* set1, const int set1_num, const int* set2, const int set2_num);

std::vector<int> linspace(int a, int b);

void newDivisionGeneration(std::vector<std::vector<int> >& D, int connect_cell_num, std::vector<int> constraint);

void switchCellInProblem(GraphProblem* P, int process_num = 0);

