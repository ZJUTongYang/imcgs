#include "demo_graph.h"
#include "graph_problem.h"
//
//GraphProblem demo1()
//{
//	std::vector<int> edge_0{ 0, 1 };
//	std::vector<int> edge_1{ 1, 2 };
//
//	std::vector<int> color_0{ 5, 6 };
//	std::vector<int> color_1{ 7, 8 };
//
//	int isRemoved_0 = -1;
//	int isRemoved_1 = -1;
//
//	Cell cell_0(edge_0, color_0, isRemoved_0);
//	Cell cell_1(edge_1, color_1, isRemoved_1);
//
//	Edge e0(-1, 0, 0); // The index of cell starts from 0 in C++
//	Edge e1(0, 1, 0);
//	Edge e2(-1, 1, 0);
//
//	GraphProblem newP;
//	newP.cell_.clear();
//	newP.cell_.push_back(cell_0);
//	newP.cell_.push_back(cell_1);
//
//	//newP.edge_.assign(edge.begin(), edge.end());
//	newP.edge_.clear();
//	newP.edge_.push_back(e0);
//	newP.edge_.push_back(e1);
//	newP.edge_.push_back(e2);
//
//	newP.cell_subdivision_process_ = 0;
//
//	newP.changes_ = " ";
//
//	return newP;
//}
//
//GraphProblem demo2()
//{
//	std::vector<int> edge_0{ 0, 2 };
//	std::vector<int> edge_1{ 1, 4 };
//	std::vector<int> edge_2{ 2, 4, 6, 5 };
//	std::vector<int> edge_3{ 5, 3 };
//	std::vector<int> edge_4{ 6, 7 };
//
//	std::vector<int> color_0{ 0 };
//	std::vector<int> color_1{ 1 };
//	std::vector<int> color_2{ 0, 1 };
//	std::vector<int> color_3{ 0 };
//	std::vector<int> color_4{ 1 };
//
//
//	int isRemoved_0 = -1;
//	int isRemoved_1 = -1;
//	int isRemoved_2 = -1;
//	int isRemoved_3 = -1;
//	int isRemoved_4 = -1;
//
//	Cell cell_0(edge_0, color_0, isRemoved_0);
//	Cell cell_1(edge_1, color_1, isRemoved_1);
//	Cell cell_2(edge_2, color_2, isRemoved_2);
//	Cell cell_3(edge_3, color_3, isRemoved_3);
//	Cell cell_4(edge_4, color_4, isRemoved_4);
//
//	Edge e0(-1, 0, 0); // The index of cell starts from 0 in C++
//	Edge e1(-1, 1, 0);
//	Edge e2(0, 2, 0);
//	Edge e3(-1, 3, 0);
//	Edge e4(1, 2, 0);
//	Edge e5(2, 3, 0);
//	Edge e6(2, 4, 0);
//	Edge e7(-1, 4, 0);
//
//
//	GraphProblem newP;
//	newP.cell_.clear();
//	newP.cell_.push_back(cell_0);
//	newP.cell_.push_back(cell_1);
//	newP.cell_.push_back(cell_2);
//	newP.cell_.push_back(cell_3);
//	newP.cell_.push_back(cell_4);
//
//	//newP.edge_.assign(edge.begin(), edge.end());
//	newP.edge_.clear();
//	newP.edge_.push_back(e0);
//	newP.edge_.push_back(e1);
//	newP.edge_.push_back(e2);
//	newP.edge_.push_back(e3);
//	newP.edge_.push_back(e4);
//	newP.edge_.push_back(e5);
//	newP.edge_.push_back(e6);
//	newP.edge_.push_back(e7);
//
//	newP.cell_subdivision_process_ = 0;
//
//	newP.changes_ = " ";
//
//
//	return newP;
//}