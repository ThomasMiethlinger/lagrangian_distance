// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include <list> // std::list
#include <utility> // std::pair
#include <vector> // std::vector

#pragma once

namespace util
{
    namespace general
    {
        int task_index_from_time_indices(int i, int j, int n);

        int mark_ready_task(
            std::vector<std::pair<int, int>>& total_task_vector,
            std::list<std::pair<int, int>>& total_ready_time_pairs,
            std::vector<int>& total_task_indices
        );

        std::vector<int> create_task_vector(int ntask, int world_rank, int world_size);

        std::vector<std::vector<int>> partition_total_task_vector(int ntask, int world_size);

        std::vector<std::vector<int>> partition_total_task_vector(std::vector<int> task_vector, int world_size);
    }
};
