#include "graph_problem.h"
#include <stack>
#include <vector>
#include <algorithm>
#include <iterator>
#include <string>
#include <iostream>
#include <list>
#include <Windows.h>
#include "commons.h"

// We pre-allocate some static variables, so that we needn't delete them in every loop
// So they should be of the maximum possible size. 
static int edge_num[NUM_THREAD + 1];
static int connect_edge[NUM_THREAD + 1][MAX_EDGE_NUM_PER_CELL];
static int connect_cell[NUM_THREAD + 1][MAX_EDGE_NUM_PER_CELL];
static int possible_color[NUM_THREAD + 1][MAX_COLOR_NUM_PER_CELL];
static int static_enums[NUM_THREAD + 1][2000000];
static int adj_color1[NUM_THREAD + 1][MAX_COLOR_NUM_PER_CELL];
static int adj_color2[NUM_THREAD + 1][MAX_COLOR_NUM_PER_CELL];
static int sub_inter_temp[NUM_THREAD + 1][MAX_COLOR_NUM_PER_CELL];
static int sub_inter[NUM_THREAD + 1][MAX_COLOR_NUM_PER_CELL];

void manipulateProblem(int id)
{
	thread_finish[id] = false;
	// A combination of the whole process, for multi-thread implementation
	int local_min_cost = 100;

	std::stack<GraphProblem*> local_S;
	
	while (1)
	{
		GraphProblem* oldP_end;

		if (stop_all_threads)
		{
			// We send all unsolved problem to the global stack, so that we can see how many sub-problems remain unsolved
			if (!local_S.empty())
			{
				S_mtx.lock();
				while (!local_S.empty())
				{
					GraphProblem* p = local_S.top();
					S.push(p);
					local_S.pop();
				}
				S_mtx.unlock();
			}
			break;

		}

		if (!local_S.empty())
		{
			oldP_end = local_S.top();
			local_S.pop();
			if (empty_threads && !local_S.empty())
			{
				if (S_mtx.try_lock())
				{
					double n = local_S.size();
					int num_to_be_sent = floor((double)n * NUM_THREAD / (NUM_THREAD + 1));
					for (unsigned int i = 0; i < num_to_be_sent; ++i)
					{
						if (local_S.empty())
							break;

						auto& temp = local_S.top();
						S.push(temp);
						local_S.pop();
					}
					S_mtx.unlock();
					empty_threads = false;
				}
			}
		}
		else
		{
			thread_finish[id] = true;
			if (S_mtx.try_lock())
			{
				if (!S.empty())
				{
					oldP_end = S.top();
					S.pop();

					int n = S.size();
					int num_to_be_get = floor((double)n / (NUM_THREAD + 1));
					for (unsigned int i = 0; i < num_to_be_get; ++i)
					{
						local_S.push(S.top());
						S.pop();
					}
					thread_finish[id] = false;
					S_mtx.unlock();
				}
				else
				{
					empty_threads = true;
					S_mtx.unlock();
					Sleep(1);
					continue;
				}
			}
			else
			{
				Sleep(1);
				continue;
			}
		}

		if (oldP_end->cell_subdivision_process_ >= oldP_end->num_all_cell_)
		{
			local_leaf_cases[id]++;
			if(local_leaf_cases[id] % 1000 == 0)
				std::cout << "Thread [" << id << "]: local leaf case: " << local_leaf_cases[id] << std::endl;

			int the_min_cost;
			std::pair<GraphProblem*, std::list<int*> > the_optimal_solution;
			IntersectionFreeGraphPainter(oldP_end, the_min_cost, the_optimal_solution, local_min_cost, id);

			if (local_min_cost < global_min_cost)
			{
				global_min_cost_mtx.lock();
				global_min_cost = local_min_cost;
				global_min_cost_mtx.unlock();

			}
			else if (local_min_cost > global_min_cost)
			{
				local_min_cost = global_min_cost;
				// We cannot directly clear them, because they store the pointers
				for (auto iter = local_optimal_solution_list[id].begin(); iter != local_optimal_solution_list[id].end(); ++iter)
				{
					iter->first->clearMemory();
					delete iter->first;
					for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2)
					{
						delete[] * iter2;
					}
				}
				local_optimal_solution_list[id].clear();
			}

			if (local_min_cost > the_min_cost)
			{
				for (auto iter = local_optimal_solution_list[id].begin(); iter != local_optimal_solution_list[id].end(); ++iter)
				{
					iter->first->clearMemory();
					delete iter->first;
					for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2)
					{
						delete[] * iter2;
					}
				}
				local_optimal_solution_list[id].clear();
				local_optimal_solution_list[id].emplace_back(the_optimal_solution);
				local_min_cost = the_min_cost;
			}
		}
		else
		{
			divideCell(oldP_end, local_S, id);
			local_total_cases[id]++;
		}
	}

	std::cout << "Thread: " << id << " return " << std::endl;
}

void divideCell(GraphProblem* P, std::stack<GraphProblem*>& local_S, const int id)
{
	int cell_index = P->cell_subdivision_process_;
	if (P->cell_[cell_index].isRemoved_ == 1)
	{
		P->cell_subdivision_process_ += 1;
		local_S.push(P);
		return ;
	}

	edge_num[id] = P->cell_[cell_index].num_edge_;

	memcpy(connect_edge[id], P->cell_[cell_index].connect_edge_, sizeof(int)*edge_num[id]);

	for (unsigned int i = 0; i < edge_num[id]; ++i)
	{
		connect_cell[id][i] = (P->edge_[connect_edge[id][i]*3] == cell_index) ? P->edge_[connect_edge[id][i]*3+1] : P->edge_[connect_edge[id][i]*3];
	}

	// Get possible color
	int color_num = P->cell_[cell_index].num_color_;

	// If this cell only have one possible color, we cannot divide it further
	if (color_num == 1)
	{
		P->cell_subdivision_process_ += 1;
		local_S.push(P);
		return;
	}

	memcpy(possible_color[id], P->cell_[cell_index].possible_color_, sizeof(int)*color_num);

	for(unsigned int i = 0; i < edge_num[id]; ++i)
	{
		int theOtherCell = connect_cell[id][i];
		if (theOtherCell != -1 && P->edge_[connect_edge[id][i]*3+2] == 1)
		{
			if(isEmptyColorIntersect(possible_color[id], color_num, P->cell_[theOtherCell].possible_color_, P->cell_[theOtherCell].num_color_))
			{
				P->clearMemory();
				delete P;
				return ;
			}
		}
	}

	// If there is no more than 3 edges that we can freely assigned, don't do sub-division
	if (P->isIntersectionFree(cell_index))
	{
		P->cell_subdivision_process_ += 1;
		local_S.push(P);
		return ;
	}

	// Here the paremeter is related to this cell
	std::vector<int> thisConstraint;
	for (unsigned int i = 0; i < edge_num[id]; ++i)
	{
		if (P->edge_[connect_edge[id][i]*3+2] == -1)
		{
			// Here we should not exclude the adjacent cell "-1". 
			// Just in opposite, we aim to include them as many as possible
			// We can connect the cell "-1", because leaving them to the main sub-cell is equivalent to connect them. 
			thisConstraint.push_back(i);
		}
	}

	std::vector<std::vector<int> > D;
	newDivisionGeneration(D, edge_num[id], thisConstraint);

	// Here we should not consider thisEnforce, because enforcing connection 
	// does not mean that it is enforced to be connected to the original cell 
	// (may be connected to the new sub-cell). 
	for (unsigned int i = 0; i < D.size(); ++i)
	{
		std::vector<std::vector<int> > subCellIndices;

		CellPartition(subCellIndices, &(D[i][0]), connect_cell[id], edge_num[id]);

		GraphProblem* newPtemp = new GraphProblem(*P);

		if (subCellIndices.empty())
		{
			newPtemp->cell_subdivision_process_ += 1;
			local_S.push(newPtemp);
			continue;
		}

		bool unreasonable_partition = false;

		// This records the index of newly generated sub-cell (for debug)
		std::vector<int> indices_of_new_nodes;


		for (unsigned int j = 0; j < subCellIndices.size(); ++j)
		{
			if (newPtemp->cell_[3].num_color_ == 0)
				coutProblem(newPtemp);

			if (subCellIndices[j].size() == 1)
			{
				// Since we will discard the "constraint" on
				// edges when painting, hence it is unnecessary to create
				// another problem
				unreasonable_partition = true;
				break;
			}
			else
			{

				// If the two adjancet cells cannot be connectable, then we
				// shouldn't partition this cell like this, because we only
				// need to leave them

				// We check that edges should be allowed to connected
				if (connect_cell[id][subCellIndices[j].front()] == -1 ||
					connect_cell[id][subCellIndices[j].back()] == -1)
				{
					unreasonable_partition = true;
					// because we have another partitioning that includes the adjacent cell "-1" into
					// the main sub-cell
					break;
				}

				int num_adj_color1 = newPtemp->cell_[connect_cell[id][subCellIndices[j].front()]].num_color_;

				memcpy(adj_color1[id], newPtemp->cell_[connect_cell[id][subCellIndices[j].front()]].possible_color_, sizeof(int)*num_adj_color1);

				int num_adj_color2 = newPtemp->cell_[connect_cell[id][subCellIndices[j].back()]].num_color_;

				memcpy(adj_color2[id], newPtemp->cell_[connect_cell[id][subCellIndices[j].back()]].possible_color_, sizeof(int)*num_adj_color2);

				int num_sub_inter_temp = colorIntersect(possible_color[id], color_num, adj_color1[id], num_adj_color1, sub_inter_temp[id]);
				if (num_sub_inter_temp == 0)
				{
					unreasonable_partition = true;
					break;
				}

				int num_sub_inter = colorIntersect(sub_inter_temp[id], num_sub_inter_temp, adj_color2[id], num_adj_color2, sub_inter[id]);
				if (num_sub_inter == 0)
				{
					unreasonable_partition = true;
					break;
				}

				// We need to check that the exact un-chosen adjacent cells must have a different color to the chosen one
				if (num_adj_color1 == 1)
				{
					int exact_pre_index = (subCellIndices[j][0] - 1 >= 0) ? subCellIndices[j][0] - 1 : edge_num[id] -1;
					int exact_pre = connect_cell[id][exact_pre_index];
					if (exact_pre != -1 &&
						newPtemp->cell_[exact_pre].num_color_ == 1 &&
						newPtemp->cell_[exact_pre].possible_color_[0] == adj_color1[id][0])
					{
						unreasonable_partition = true;
						break;
					}
				}
				if (num_adj_color2 == 1)
				{
					int exact_post_index = (subCellIndices[j].back() + 1 < edge_num[id]) ? subCellIndices[j].back() + 1 : 0;
					int exact_post = connect_cell[id][exact_post_index];
					if (exact_post != -1 &&
						newPtemp->cell_[exact_post].num_color_ == 1 &&
						newPtemp->cell_[exact_post].possible_color_[0] == adj_color2[id][0])
					{
						unreasonable_partition = true;
						break;
					}
				}

				// This means that the two adjacent cells must be connected
				int sizeC = newPtemp->num_all_cell_;
				int sizeE = newPtemp->num_all_edge_;

				// We create the new edge
				{
					newPtemp->edge_[newPtemp->num_all_edge_ * 3] = cell_index;
					newPtemp->edge_[newPtemp->num_all_edge_ * 3 + 1] = sizeC;
					newPtemp->edge_[newPtemp->num_all_edge_ * 3 + 2] = -1;

					newPtemp->num_all_edge_++;
				}

				// If the index of the adjacent cell is smaller
				// than this one, that mean it has been sub - divided and we only
				// need to colour it(remove its other possible colour).
				// Now that sub_inter is the intersection of colours, we can
				// directly use this
				if (connect_cell[id][subCellIndices[j].front()] < cell_index)
				{
					newPtemp->cell_[connect_cell[id][subCellIndices[j].front()]].assignColor(sub_inter[id], num_sub_inter);
				}
				if (connect_cell[id][subCellIndices[j].back()] < cell_index)
				{
					newPtemp->cell_[connect_cell[id][subCellIndices[j].back()]].assignColor(sub_inter[id], num_sub_inter);
				}

				// Here there may be a bug, because the sub - cell
				// may be further sub - divided into sub - cells, in which case the
				// newest sub - cell should have the same possible colour as the
				// original cell but not the sub - cell
				// We create the new cell

				int* new_connect_edge = new int[MAX_EDGE_NUM_PER_CELL];
				for (unsigned int k = 0; k < subCellIndices[j].size(); ++k)
				{
					new_connect_edge[k] = connect_edge[id][subCellIndices[j][k]];
				}
				new_connect_edge[subCellIndices[j].size()] = sizeE;

				newPtemp->insertCell(Cell(new_connect_edge, (subCellIndices[j].size() + 1), sub_inter[id], num_sub_inter, -1));

				delete[] new_connect_edge;

				newPtemp->cell_subdivision_process_ += 1;

				indices_of_new_nodes.push_back(sizeC);

				// We update the elements in the edgelist
				for(unsigned int k = 0; k < subCellIndices[j].size(); ++k)
				{
					int e = connect_edge[id][subCellIndices[j][k]];
					if (newPtemp->edge_[e*3] == cell_index)
						newPtemp->edge_[e*3] = sizeC;
					else
						newPtemp->edge_[e*3+1] = sizeC;

					if (newPtemp->edge_[e*3] > newPtemp->edge_[e*3+1])
					{
						int temp = newPtemp->edge_[e*3];
						newPtemp->edge_[e*3] = newPtemp->edge_[e*3+1];
						newPtemp->edge_[e*3+1] = temp;
					}
				}

				newPtemp->edge_[connect_edge[id][subCellIndices[j].front()]*3+2] = 1;
				newPtemp->edge_[connect_edge[id][subCellIndices[j].back()]*3+2] = 1;

				// We replace all those edge indices temporarily
				for (unsigned int k = 0; k < subCellIndices[j].size(); ++k)
				{
					newPtemp->cell_[cell_index].connect_edge_[subCellIndices[j][k]] = sizeE;
				}

				newPtemp->changes_ += std::to_string(cell_index);
				newPtemp->changes_ += "->";
				newPtemp->changes_ += std::to_string(newPtemp->num_all_cell_ - 1);
				newPtemp->changes_ += ",";

			} // for each valid subCellIndices
		} // for each subCellIndices

		if (unreasonable_partition)
		{
			newPtemp->clearMemory();
			delete newPtemp;
			continue;
		}

		// We update the edgelist of the old cell
		std::vector<int> temp;
		temp.resize(newPtemp->cell_[cell_index].num_edge_, 0);
		for (unsigned int k = 0; k < newPtemp->cell_[cell_index].num_edge_; ++k)
		{
			temp[k] = newPtemp->cell_[cell_index].connect_edge_[k];
		}

		while (1)
		{
			if (temp.empty() || temp.front() != temp.back())
			{
				break;
			}
			temp.erase(temp.end() - 1);
		}
		for (int j = temp.size() - 2; j >= 0; --j)
		{
			if (temp[j] == temp[j + 1])
				temp.erase(temp.begin() + j);
		}

		newPtemp->cell_[cell_index].assignEdge(&(temp[0]), temp.size());

		// We enum all colours of this cell
		const long long int bit_1 = 1;
		long long int bit_inter = 0;
		for (unsigned int k = 0; k < color_num; ++k)
		{
			bit_inter |= bit_1 << possible_color[id][k];
		}

		// The possible colour cannot be the same as the adjacent ones 
		for (unsigned int j = 0; j < indices_of_new_nodes.size(); ++j)
		{
			int c = indices_of_new_nodes[j];
			if (newPtemp->cell_[indices_of_new_nodes[j]].num_color_ == 1)
			{
				bit_inter ^= bit_1 << newPtemp->cell_[c].possible_color_[0];
			}
		}

		// If there is no possible colour to use
		if (bit_inter == 0)
		{
			newPtemp->clearMemory();
			delete newPtemp;
			continue;
		}

		int* inter = new int[MAX_COLOR_NUM_PER_CELL];
		int inter_num = 0;
		for (unsigned int k = 0; k < color_num; ++k)
		{
			if (bit_inter & bit_1 << possible_color[id][k])
			{
				inter[inter_num] = possible_color[id][k];
				inter_num++;
			}
		}

		newPtemp->cell_subdivision_process_ += 1;

		newPtemp->cell_[cell_index].assignColor(inter, inter_num);

		local_S.push(newPtemp);
		delete[] inter;

	} // for each element in D

	P->clearMemory();
	delete P;
}



void newCollectIntersectionFreeCells(GraphProblem* P, int* inter_free, std::vector<std::vector<int> >& adj_cells_list)
{

	// For easy coding, we assume that cell 1 cannot be an internal cell
	// (Or else, we need to find the index of the first cell to be coloured)
	inter_free[0] = -1;
	std::vector<int> temp;
	adj_cells_list.resize(P->num_all_cell_, temp);

	for (unsigned int i = 1; i < P->num_all_cell_; ++i)
	{
		// If the cell has only one possible colour, it need not be an internal cell
		// (Should it be an internal cell, it would prevent its adjacent cells to be an internal cell. )
		if (inter_free[i] == -1)
		{
			continue;
		}
		else if (inter_free[i] == 1)
		{
			// This is impossible
		}
		else // if inter_free[i] == 0
		{
			if (P->cell_[i].num_color_ == 1)
			{
				inter_free[i] = -1;
				continue;
			}
			// we check whether all its connecable adjacent cells have been assigned 
			std::vector<int> A;
			P->adjConnectableCells(A, i);
			bool non_free = false;

			if (A.size() >= 4)
				non_free = true;

			for (auto iter = A.begin(); iter != A.end(); ++iter)
			{
				if (inter_free[*iter] > 0)
				{
					non_free = true;
					break;
				}
			}
			if (non_free)
			{
				inter_free[i] = -1;
			}
			else
			{
				inter_free[i] = 1;
				for (auto iter = A.begin(); iter != A.end(); ++iter)
				{
					inter_free[*iter] = -1;
				}
				adj_cells_list[i].swap(A);
			}
		}
	}
}

struct Edge {
	int c1_;
	int c2_;
	int constraint_;
	Edge()
	{
		constraint_ = 0;
	}
	Edge(int c1, int c2, int cons)
	{
		c1_ = c1;
		c2_ = c2;
		constraint_ = cons;
	}
};

void IntersectionFreeGraphPainter(GraphProblem* P, int& min_cost, std::pair<GraphProblem*, std::list<int*> >& optimal_solution, const int local_min_cost, const int id)
{
	// In this function, we try to assign colours based on the edge
	const long long int bit_1 = 1;

	// We firstly find all intersection-free cells
	// For these cells, they will have the exactly ONE solution, so we need not enumerate them when creating the solution list

	int* inter_free = new int[P->num_all_cell_];
	memset(inter_free, 0, sizeof(int)*P->num_all_cell_);

	std::vector<std::vector<int> > adj_cells_list;

	newCollectIntersectionFreeCells(P, inter_free, adj_cells_list);

	// We filter all edges. We preseve the useful edges, and order them. 
	// We will use them when creating the painting sequences, instead of checking afterwards
	std::vector<Edge> X;
	{
		Edge temp;
		X.resize(P->num_all_edge_, temp);
		for (unsigned int i = 0; i < P->num_all_edge_; ++i)
		{
			X[i].c1_ = P->edge_[i * 3];
			X[i].c2_ = P->edge_[i * 3 + 1];
			X[i].constraint_ = P->edge_[i * 3 + 2];
		}
	}

	X.erase(std::remove_if(X.begin(), X.end(), [&](const Edge& a)
		{return a.c1_ == -1 || a.c2_ == -1 || a.constraint_ == 0 || inter_free[a.c1_] == 1 || inter_free[a.c2_] == 1; }
		), X.end());

	std::sort(X.begin(), X.end(), [](const Edge& a, const Edge& b) 
		{return a.c1_ < b.c1_ || (a.c1_ == b.c1_ && a.c2_ < b.c2_); }
		);

	int num_of_cell = P->num_all_cell_;

	int max_num_of_enums = 1;
	for (unsigned int i = 0; i < num_of_cell; ++i)
	{
		if(inter_free[i] != 1)
			max_num_of_enums *= P->cell_[i].num_color_;
	}

	// We enumerate all possible combinations of cell colours
	// Even if a cell is isolated, we assign a colour to it
	// If the cell is internal, then we do not assign it (let it be -1)
	// We should filter the constraints before the internal cells are assigned, 
	// or else we may remove some coincidently constraint-violated internal cells
	int count = 0;

	bool use_static_memory = false;
	int* enums;
	if (max_num_of_enums < 2000000)
	{
		enums = static_enums[id];
		use_static_memory = true;
	}
	else
	{
		enums = new int[max_num_of_enums*num_of_cell];
	}

	for (unsigned int i = 0; i < P->cell_[0].num_color_; ++i)
	{
		enums[i*num_of_cell] = P->cell_[0].possible_color_[i];
	}
	count = P->cell_[0].num_color_;

	//We set the possible colours for other cells
	for (unsigned int cell_index = 1; cell_index < P->num_all_cell_; ++cell_index)
	{
		if (inter_free[cell_index] == 1)
		{
			for (unsigned int i = 0; i < count; ++i)
			{
				enums[i*num_of_cell + cell_index] = -1;
			}
		}
		else
		{
			int first_color = P->cell_[cell_index].possible_color_[0];
			for (unsigned int i = 0; i < count; ++i)
			{
				enums[i*num_of_cell + cell_index] = first_color;
			}

			if (P->cell_[cell_index].num_color_ > 1)
			{
				for (unsigned int color_index = 1; color_index < P->cell_[cell_index].num_color_; ++color_index)
				{
					memcpy(&(enums[color_index*count*num_of_cell]), enums, sizeof(int)*count*num_of_cell);
					for (unsigned int i = 0; i < count; ++i)
					{
						enums[(color_index*count + i)*num_of_cell + cell_index] = P->cell_[cell_index].possible_color_[color_index];
					}
				}
				count *= P->cell_[cell_index].num_color_;
			}
		}

		int count_minus_one = count - 1;
		while (!X.empty() && X[0].c2_ <= cell_index)
		{
			auto iter = X.begin();
			int i = 0;
			while (i < count_minus_one - 1)
			{
				if ( (iter->constraint_ == 1 && (bit_1 << enums[i*num_of_cell + iter->c1_] ^ bit_1 << enums[i*num_of_cell + iter->c2_])) || 
					 (iter->constraint_ == -1 && (bit_1 << enums[i*num_of_cell + iter->c1_] & bit_1 << enums[i*num_of_cell + iter->c2_])) )
				{
					// The colour enumeration case i is not valid, we find the last valid item and move it here
					while (count_minus_one > i + 1 &&	
						((iter->constraint_ == 1 && (enums[count_minus_one*num_of_cell + iter->c1_] != enums[count_minus_one*num_of_cell + iter->c2_])) ||	
						(iter->constraint_ == -1 && (enums[count_minus_one*num_of_cell + iter->c1_] == enums[count_minus_one*num_of_cell + iter->c2_]))))
					{
						--count_minus_one;
					}
					
					if (count_minus_one > i + 1)
					{
						memcpy(&(enums[i*num_of_cell]), &(enums[count_minus_one*num_of_cell]), sizeof(int)*num_of_cell);
						--count_minus_one;
					}
					else
					{
						break;
					}
				}
				++i;
			}

			X.erase(X.begin());
		}//while
		count = count_minus_one + 1;
	}

	delete[] inter_free;
	inter_free = nullptr;

	int current_best_cost = local_min_cost;
	for (unsigned int i = 0; i < count; ++i)
	{
		// Even if there are enums that violates the constraint, the number of them is small (should be < 4). So we directly enumerate them. 

		int c = P->graphPainter(&(enums[i*num_of_cell]), adj_cells_list, local_min_cost);
		optimal_solution.first = P;
		if (c < current_best_cost)
		{
			current_best_cost = c;
			optimal_solution.second.clear();
			int* result = new int[num_of_cell];
			memcpy(result, &(enums[i*num_of_cell]), sizeof(int)*num_of_cell);
			optimal_solution.second.emplace_back(result);
		}
		else if (c == current_best_cost)
		{
			int* result = new int[num_of_cell];
			memcpy(result, &(enums[i*num_of_cell]), sizeof(int)*num_of_cell);
			optimal_solution.second.emplace_back(result);
		}
	}

	if (optimal_solution.second.empty())
	{
		P->clearMemory();
		delete P;
		optimal_solution.first = nullptr;
	}

	// Since the min_cost is the local min cost at the beginning, so we might return an empty structure.
	// So we must check that optimal_solution.first != nullptr before storing the solution
	min_cost = current_best_cost;

	if(!use_static_memory)
		delete[] enums;
}

int GraphProblem::graphPainter(int* S, const std::vector<std::vector<int> >& adj_cells_list, const int local_min_cost, bool check)
{
	// If there is some colour -1, that means it is an internal cell. We just assign its colour
	for (unsigned int i = 0; i < num_all_cell_; ++i)
	{
		if (S[i] == -1)
		{
			if (adj_cells_list[i].empty())
				S[i] = cell_[i].possible_color_[0];
			else
				S[i] = bestColor(&(adj_cells_list[i][0]), adj_cells_list[i].size(), &(cell_[i].possible_color_[0]), cell_[i].num_color_, S);
		}
	}
	
	// We should check that this painting does not use more than min-number of colours. 
	// Or else we can directly remove it
	int num_of_used_colours = cal_num_of_used_colours(S, num_all_cell_);
	if (num_of_used_colours > local_min_cost)
	{
		return 101;
	}

	// Disconnectable regions will be assigned different colours
	int n = num_all_cell_;
	int* painted = new int[n];
	memset(painted, -1, sizeof(int)*n);

	// For the edges which connect two same-colour cells, we preserve it, or else we remove it. 
	//std::vector<Edge> X(edge_.begin(), edge_.end());
	std::vector<Edge> X;
	{
		Edge temp;
		X.resize(num_all_edge_, temp);
		for (unsigned int i = 0; i < num_all_edge_; ++i)
		{
			X[i].c1_ = edge_[i * 3];
			X[i].c2_ = edge_[i * 3 + 1];
			X[i].constraint_ = edge_[i * 3 + 2];
		}
	}

	// Here we cannot directly manipulate S, because cells with the same colour might not be continuously trackable
	if (check)
		coutProblem(this);

	if (check)
		coutAVector(S, num_all_cell_);

	if (check)
	{
		for (auto iter = X.begin(); iter != X.end(); ++iter)
			std::cout << "[" << iter->c1_ << ", " << iter->c2_ << ", " << iter->constraint_ << "]" << std::endl;
	}

	for (auto iter = X.begin(); iter != X.end(); ++iter)
	{
		if (iter->c1_ == -1 || iter->c2_ == -1 || S[iter->c1_] != S[iter->c2_] || iter->constraint_ == -1)
			continue;

		if (painted[iter->c1_] == -1 && painted[iter->c2_] == -1)
		{
			painted[iter->c1_] = iter - X.begin();
			painted[iter->c2_] = iter - X.begin();
		}
		else if (painted[iter->c1_] == -1 && painted[iter->c2_] != -1)
		{
			painted[iter->c1_] = painted[iter->c2_];
		}
		else if (painted[iter->c1_] != -1 && painted[iter->c2_] == -1)
		{
			painted[iter->c2_] = painted[iter->c1_];
		}
		else if (painted[iter->c1_] != painted[iter->c2_])
		{
			for (unsigned int i = 0; i < n; ++i)
			{
				painted[i] = (painted[i] == painted[iter->c2_]) ? painted[iter->c1_] : painted[i];
			}
		}

		if (check)
		{
			std::cout << "[" << iter->c1_ << ", " << iter->c2_ << ", " << iter->constraint_ << "]" << std::endl;
			std::cout << "check painted: ";
			coutAVector(painted, num_all_cell_);
		}
	}

	if (check)
	{
		std::cout << "we check the first time" << std::endl;
		coutAVector(painted, n);
	}

	// Since we removed some cell, there might be some uncoloured cells. We label them as new indices
	int max_element = -1;
	for (unsigned int i = 0; i < n; ++i)
	{
		max_element = (max_element < painted[i]) ? painted[i] : max_element;
	}
	max_element += 1;
	
	if(check)
		std::cout << "max_element = " << max_element << std::endl;
	
	for (unsigned int i = 0; i < n; ++i)
	{
		if (painted[i] == -1)
			painted[i] = i + max_element;
	}

	if(check)
		coutAVector(painted, n);

	// We find the number of distinct colours, which is the cost
	// Here we assume that the indice of colours are not bigger than 64
	long long int bit_colour = 0;
	const long long int bit_1 = 1;
	int cost = 0;
	for (unsigned int i = 0; i < n; ++i)
	{
		if (!(bit_colour & (bit_1 << painted[i])))
		{
			bit_colour |= bit_1 << painted[i];
			cost += 1;
		}
	}

	delete[] painted;
	return cost;
}


bool GraphProblem::isIntersectionFree(int cell_index)
{
	int count = 0;

	for (unsigned int j = 0; j < cell_[cell_index].num_edge_; ++j)
	{
		if (edge_[cell_[cell_index].connect_edge_[j] * 3] == -1 || edge_[cell_[cell_index].connect_edge_[j] * 3 + 1] == -1)
		{
			continue;
		}
		int connect_cell_j = (edge_[cell_[cell_index].connect_edge_[j] * 3] == cell_index) ? edge_[cell_[cell_index].connect_edge_[j] * 3 + 1] : edge_[cell_[cell_index].connect_edge_[j] * 3];

		if (isEmptyColorIntersect(cell_[cell_index].possible_color_, cell_[cell_index].num_color_, cell_[connect_cell_j].possible_color_, cell_[connect_cell_j].num_color_))
		{
			continue;
		}
		if (++count > 3)
		{
			return false;
		}
	}

	return true;
}

void GraphProblem::adjConnectableCells(std::vector<int>& A, int cell_index)
{
	A.reserve(cell_[cell_index].num_edge_);

	for (unsigned int j = 0; j < cell_[cell_index].num_edge_; ++j)
	{
		if (edge_[cell_[cell_index].connect_edge_[j]*3] == -1 || edge_[cell_[cell_index].connect_edge_[j]*3+1] == -1)
		{
			continue;
		}
		int connect_cell_j = (edge_[cell_[cell_index].connect_edge_[j]*3] == cell_index) ? edge_[cell_[cell_index].connect_edge_[j]*3+1] : edge_[cell_[cell_index].connect_edge_[j]*3];

		if (isEmptyColorIntersect(cell_[cell_index].possible_color_, cell_[cell_index].num_color_, cell_[connect_cell_j].possible_color_, cell_[connect_cell_j].num_color_))
		{
			continue;
		}

		A.push_back(connect_cell_j);
	}
}



