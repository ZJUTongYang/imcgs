//#include <vld.h>
#include <stdio.h>
#include "graph_problem.h"
#include <string.h>
#include <iostream>
#include <stack>
#include "graph_problem.h"
#include <thread>
#include <vector>
#include <list>
#include <algorithm>
#include <Windows.h>
#include <fstream>
#include <sstream>
#include <cassert>
#include <mutex>
#include <atomic>
#include "demo_graph.h"
#include "commons.h"

// We define some global element

volatile bool stop_all_threads;
std::atomic<bool> empty_threads;

std::mutex S_mtx;
std::stack<GraphProblem*> S; // The storage of all unsolved (unbranched or unpainted) problem

// Local Variables
long long int local_total_cases[NUM_THREAD];
long long int local_leaf_cases[NUM_THREAD];
std::list< std::pair<GraphProblem*, std::list<int*> > > local_optimal_solution_list[NUM_THREAD];

std::atomic<bool> thread_finish[NUM_THREAD];

std::mutex global_min_cost_mtx;
int global_min_cost = 100;

int main() {

	printf("We create the problem.\n");

	long long int the_total_cases = 0;
	long long int the_leaf_cases = 0;
	long long int the_optimal_leaf_cases = 0;
	int local_min_cost = 100;
	std::list< std::pair<GraphProblem*, std::list<int*> > > the_optimal_solution_list;

	GraphProblem* newP = demo_read_from_file("problem_input_cardoor_reduced.txt");

	printf("We show the original problem before cell switching: \n");
	coutProblem(newP);

	switchCellInProblem(newP);

	std::cout << "Check the modified graph with the : " << sizeof(newP) << std::endl;
	coutProblem(newP);

	GraphProblem* oldP = new GraphProblem;
	oldP->assign(*newP);

	newP->clearMemory();
	delete newP;

	S.push(oldP);

	stop_all_threads = false;
	empty_threads = false;

	for (unsigned int i = 0; i < NUM_THREAD; ++i)
	{
		thread_finish[i] = true;
	}

	long start, finish;
	long painting_time = 0, branching_time = 0;

	std::thread threads[NUM_THREAD];
	for (unsigned int i = 0; i < NUM_THREAD; ++i)
	{
		threads[i] = std::thread(manipulateProblem, i);
		threads[i].detach();
	}

	std::stack<GraphProblem*> local_S;
	start = GetTickCount();

	while (1)
	{
		GraphProblem* oldP_end;

		// If all threads have finished the task, we return. 
		bool all_finished = local_S.empty();
		
		if (all_finished)
		{
			S_mtx.lock();
			all_finished &= S.empty();
			S_mtx.unlock();
		}

		if (all_finished)
		{
			for (unsigned int i = 0; i < NUM_THREAD; ++i)
			{
				all_finished &= thread_finish[i];
			}
		}
		if (all_finished)
			break;

		//////////////////////////////////////////////////////////////////
		// We pick up an unsolved problem. 
		// If the local queue have been empty, then we pick up one from the global queue (with thread safety)
		if (!local_S.empty())
		{
			oldP_end = local_S.top();
			local_S.pop();
			if (empty_threads && !local_S.empty())
			{
				// Calculate how many cases that we should move from local queue to global queue
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

		// Start the enumeration: sub-dividing a cell, or painting the graph
		if (oldP_end->cell_subdivision_process_ >= oldP_end->num_all_cell_)
		{

			if(the_leaf_cases % 1000 == 0)
				std::cout << "local leaf cases = " << the_leaf_cases << std::endl;
	
			the_leaf_cases++;

			int the_min_cost;
			std::pair<GraphProblem*, std::list<int*> > the_optimal_solution;
			IntersectionFreeGraphPainter(oldP_end, the_min_cost, the_optimal_solution, local_min_cost);

			
			if (local_min_cost < global_min_cost)
			{
				global_min_cost_mtx.lock();
				global_min_cost = local_min_cost;
				global_min_cost_mtx.unlock();
			}
			else if (local_min_cost > global_min_cost)
			{
				local_min_cost = global_min_cost;
				for (auto iter = the_optimal_solution_list.begin(); iter != the_optimal_solution_list.end(); ++iter)
				{
					iter->first->clearMemory();
					delete iter->first;
					for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2)
					{
						delete[] * iter2;
					}
				}
				the_optimal_solution_list.clear();
			}

			if (local_min_cost > the_min_cost)
			{
				the_optimal_leaf_cases = 1;
				
				for (auto iter = the_optimal_solution_list.begin(); iter != the_optimal_solution_list.end(); ++iter)
				{
					iter->first->clearMemory();
					delete iter->first;
					for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2)
					{
						delete[] * iter2;
					}
				}
				the_optimal_solution_list.clear();
				the_optimal_solution_list.emplace_back(the_optimal_solution);
				local_min_cost = the_min_cost;
			}
			else if (local_min_cost == the_min_cost && the_optimal_solution.first != nullptr)
			{
				the_optimal_leaf_cases++;
				// Here for the memory efficiency, we only store one optimal solution

				//all_optimal_solutions.emplace_back(the_optimal_solution);
			}
		}
		else
		{
			divideCell(oldP_end, local_S);
			the_total_cases++;
		}
	}


	finish = GetTickCount();

	stop_all_threads = true;
	std::cout << "We wait for all threads return" << std::endl;
	Sleep(2000);

	std::cout << "We join the optimal solutions collected by different thread." << std::endl;
	for (unsigned int i = 0; i < NUM_THREAD; ++i)
	{
		the_optimal_solution_list.splice(the_optimal_solution_list.end(), local_optimal_solution_list[i]);
	}

	//// Here we see how many unsolved problems are still in the stack
	//std::cout << "We check all remaining unsolved problems" << std::endl;
	//if (!local_S.empty())
	//{
	//	S_mtx.lock();
	//	while (!local_S.empty())
	//	{
	//		GraphProblem* p = local_S.top();
	//		S.push(p);
	//		local_S.pop();
	//	}
	//	S_mtx.unlock();
	//}
	//std::cout << "number of unsolved problems: " << S.size() << std::endl;
	//std::vector<int> steps;
	//while (!S.empty())
	//{
	//	steps.emplace_back(S.top()->cell_subdivision_process_);
	//	S.pop();
	//}
	//std::sort(steps.begin(), steps.end());
	//coutAVector(steps);

	std::cout << "We print one of the optimal solution: " << std::endl;
	if (!the_optimal_solution_list.empty())
	{
		auto check = the_optimal_solution_list.front();
		auto optimal_colouring = check.second.front();
		coutProblem(check.first);
		coutAVector(optimal_colouring, check.first->num_all_cell_);

		std::vector<std::vector<int> > temp;
	}
	else
	{
		std::cout << "There is no optimal solution to be shown" << std::endl;
	}

	std::cout << "Total time: " << finish - start << "ms " << std::endl;

	long long int total_cases = the_total_cases;
	for (unsigned int i = 0; i < NUM_THREAD; ++i)
	{
		total_cases += local_total_cases[i];
	}
	std::cout << "Total Branched Cases: " << total_cases << std::endl;
	std::cout << "check cases in each thread: main " << the_total_cases << ", [";
	for (unsigned int i = 0; i < NUM_THREAD; ++i)
	{
		std::cout << local_total_cases[i] << ", ";
	}
	std::cout << "]" << std::endl;


	long long int leaf_cases = the_leaf_cases;
	for (unsigned int i = 0; i < NUM_THREAD; ++i)
	{
		leaf_cases += local_leaf_cases[i];
	}
	std::cout << "The number of cases that are still valid after cell sub-divisions (i.e., the number of graph paintings): " << leaf_cases << std::endl;

	std::cout << "Current optimal cost (the number of colours): " << global_min_cost << std::endl;

	// We release all optimal solutions
	for (auto iter = the_optimal_solution_list.begin(); iter != the_optimal_solution_list.end(); ++iter)
	{
		iter->first->clearMemory();
		delete iter->first;
		for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2)
		{
			delete *iter2;
		}
	}


	return 0;
}

