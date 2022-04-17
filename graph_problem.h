#pragma once

#include <vector>
#include <stack>
#include <list>
#include <mutex>
#include <iostream>
#include <atomic>

// Pre-defined variables for efficient pre-allocation of memory
#define MAX_COLOR_NUM_PER_CELL 5
#define MAX_EDGE_NUM_PER_CELL 23
#define MAX_EDGE_NUM 80

#define NUM_THREAD 11

struct GraphProblem;

extern volatile bool stop_all_threads;
extern std::atomic<bool> empty_threads;

extern std::mutex S_mtx;
extern std::stack<GraphProblem*> S;

// Local Variables
extern long long int local_total_cases[NUM_THREAD];
extern long long int local_leaf_cases[NUM_THREAD];
extern std::list< std::pair<GraphProblem*, std::list<int*> > > local_optimal_solution_list[NUM_THREAD];

extern std::atomic<bool> thread_finish[NUM_THREAD];

extern std::mutex global_min_cost_mtx;
extern int global_min_cost;

struct Cell {

	int connect_edge_[MAX_EDGE_NUM_PER_CELL];
	int possible_color_[MAX_COLOR_NUM_PER_CELL];

	int num_edge_;
	int num_color_;

	int isRemoved_; 
	Cell()
	{
		
	}
	Cell(const int* edge_index, const int edge_num, const int* colors, const int color_num, int isRemoved)
	{
		memcpy(connect_edge_, edge_index, sizeof(int)*edge_num);

		memcpy(possible_color_, colors, sizeof(int)*color_num);
		num_edge_ = edge_num;
		num_color_ = color_num;
		isRemoved_ = isRemoved;
	}
	void assign(const Cell& a)
	{
		memcpy(connect_edge_, a.connect_edge_, sizeof(int)*a.num_edge_);

		memcpy(possible_color_, a.possible_color_, sizeof(int)*a.num_color_);

		num_edge_ = a.num_edge_;
		num_color_ = a.num_color_;
		isRemoved_ = a.isRemoved_;
	}
	void assignColor(const int* new_color, const int n)
	{
		memcpy(possible_color_, new_color, sizeof(int)*n);
		num_color_ = n;
	}
	void assignEdge(const int* new_edge, const int n)
	{
		memcpy(connect_edge_, new_edge, sizeof(int)*n);
		num_edge_ = n;
	}

	Cell(const Cell& a)
	{
		memcpy(connect_edge_, a.connect_edge_, sizeof(int)*a.num_edge_);
	
		memcpy(possible_color_, a.possible_color_, sizeof(int)*a.num_color_);

		num_edge_ = a.num_edge_;
		num_color_ = a.num_color_;
		isRemoved_ = a.isRemoved_;
	}
	Cell& operator=(const Cell&) = delete;

};


struct GraphProblem {
	// In this new structure 
	// Given an n-edge cell, the maximum number of sub-cells is n-2. 
	// So given the graph, we can pre-calculate the maximum possible memory. 

	Cell* cell_;
	int edge_[MAX_EDGE_NUM * 3]; // For each row, three elements: [adjacent cell 1, adjacent cell 2, constraint (-1 for disconnectable, 1 for connectable)]

	int num_all_cell_;
	int num_all_edge_;

	int cell_subdivision_process_; // The index of the next cell to be manipulated. If it is >= cell_num that means we have finished subdividing all cells
	std::string changes_;

	GraphProblem():cell_(nullptr)//, edge_(nullptr)
	{
		memset(edge_, 0, sizeof(int)*MAX_EDGE_NUM * 3);
	}
	void assign(const GraphProblem& a)
	{
		if (cell_ != nullptr)
		{
			delete[] cell_;
			cell_ = nullptr;
		}
		cell_ = new Cell[a.num_all_cell_];
		memcpy(cell_, a.cell_, sizeof(Cell)*a.num_all_cell_);

		memcpy(edge_, a.edge_, sizeof(int)*a.num_all_edge_*3);

		num_all_cell_ = a.num_all_cell_;
		num_all_edge_ = a.num_all_edge_;

		cell_subdivision_process_ = a.cell_subdivision_process_;
		changes_ = a.changes_;
	}
	void insertCell(const Cell& c)
	{
		Cell* temp = cell_;
		cell_ = new Cell[num_all_cell_ + 1];
		memcpy(cell_, temp, sizeof(Cell)*num_all_cell_);
		cell_[num_all_cell_].assign(c);
		delete[] temp;
		num_all_cell_ += 1;
	}
	GraphProblem(const GraphProblem& a)
	{
		cell_ = new Cell[a.num_all_cell_];

		memcpy(cell_, a.cell_, sizeof(Cell)*a.num_all_cell_);
		memcpy(edge_, a.edge_, sizeof(int)*a.num_all_edge_ * 3);

		num_all_cell_ = a.num_all_cell_;
		num_all_edge_ = a.num_all_edge_;
		
		cell_subdivision_process_ = a.cell_subdivision_process_;
		changes_ = a.changes_;
	}
	inline void clearMemory()
	{
		if (cell_ != nullptr)
		{
			delete[] cell_;
			cell_ = nullptr;
		}
	}

	~GraphProblem()
	{
		if(cell_ != nullptr)
		{
			delete[] cell_;
			cell_ = nullptr;
		}
	}
	GraphProblem& operator=(const GraphProblem&) = delete;

	bool isIntersectionFree(int cell_index);
	void adjConnectableCells(std::vector<int>& A, int cell_index);

	int graphPainter(int* S, const std::vector<std::vector<int> >& adj_cells_list, const int local_min_cost, bool check = false);
};

void newCollectIntersectionFreeCells(GraphProblem* P, int* inter_free, std::vector<std::vector<int> >& adj_cells_list);

void manipulateProblem(int id);

void divideCell(GraphProblem* oldP_end, std::stack<GraphProblem*>& local_S, const int id = NUM_THREAD); // default parameter = main parameter = the last one

void IntersectionFreeGraphPainter(GraphProblem* P, int& min_cost, std::pair<GraphProblem*, std::list<int*> >& optimal_solution, const int local_min_cost, const int id = NUM_THREAD); // default parameter = main parameter = the last one
