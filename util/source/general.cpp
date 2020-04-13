// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include "../include/general.hpp"

#include <algorithm> // std::sort
#include <list> // std::list
#include <utility> // std::pair
#include <vector> // std::vector

using std::pair;
using std::list;
using std::vector;

typedef vector<int> VI;
typedef vector<VI> VVI;
typedef pair<int, int> PII;

int util::general::task_index_from_time_indices(int i, int j, int n)
{
    return - 1 + i * (n - 2) - (i * (i - 1)) / 2 + j;
}

bool compare_time_pairs(PII pa, PII pb)
{
    if(pa.first < pb.first)
    {
        return true;
    }
    else if(pa.first == pb.first)
    {
        if(pa.second <= pb.second)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

int util::general::mark_ready_task(vector<PII>& total_task_vector, list<PII>& total_ready_time_pairs, VI& total_task_indices)
{
    int ntask = total_task_vector.size();

    std::sort(total_task_vector.begin(), total_task_vector.end(), compare_time_pairs);
    total_ready_time_pairs.sort(compare_time_pairs);

    for(int i = 0; i < ntask; i++)
        total_task_indices[i] = i;

    int kmin = 0;
    for(auto it = total_ready_time_pairs.begin(); it != total_ready_time_pairs.end(); it++)
    {
        /*pair<pair<int, int>, double> p = *it;
        int ti = p.first.first;
        int tj = p.first.second;
        int i = (ti - tmin) / tstep;
        int j = (tj - tmin) / tstep;*/
        //total_task_indices[task_index_from_time_indices(i, j, nsteps)] = -1;
        for(int k = kmin; k < total_task_vector.size(); k++)
        {
            if((*it) == total_task_vector[k])
            {
                total_task_indices[k] = -1;
                kmin = k + 1;
            }
        }
    }
    int ctr_ready_results = 0;
    for(int i = 0; i < ntask; i++)
    {
        if(total_task_indices[i] == -1)
            ctr_ready_results++;
    }

    return ctr_ready_results;
}

VI util::general::create_task_vector(int ntask, int world_rank, int world_size)
{
    VI task_vector;
    int ntask_per_rank_min = ntask / world_size;
    int ntask_per_rank_max = ntask % world_size == 0 ? ntask_per_rank_min : ntask_per_rank_min + 1;
    int min_rank = ntask % world_size;
    int ntask_my_rank = world_rank < min_rank ? ntask_per_rank_max : ntask_per_rank_min;
    task_vector.resize(ntask_my_rank);

    int itaskmin, itaskmax;
    if(world_rank < min_rank)
    {
        itaskmin = world_rank * ntask_per_rank_max;
        itaskmax = itaskmin + ntask_per_rank_max - 1;
    }
    else
    {
        itaskmin = min_rank * ntask_per_rank_max + (world_rank - min_rank) * ntask_per_rank_min;
        itaskmax = itaskmin + ntask_per_rank_min - 1;
    }

    for(int i = itaskmin, j = 0; i <= itaskmax; i++, j++)
        task_vector[j] = i;

    return task_vector;
}

VVI util::general::partition_total_task_vector(int ntask, int world_size)
{
    VVI task_vector_partitioned(world_size);
    int ntask_per_rank_min = ntask / world_size;
    int ntask_per_rank_max = ntask % world_size == 0 ? ntask_per_rank_min : ntask_per_rank_min + 1;
    int min_rank = ntask % world_size;

    for(int r = 0; r < world_size; r++)
    {
        int ntask_my_rank = r < min_rank ? ntask_per_rank_max : ntask_per_rank_min;
        task_vector_partitioned[r].resize(ntask_my_rank);
        int itaskmin, itaskmax;

        if(r < min_rank)
        {
            itaskmin = r * ntask_per_rank_max;
            itaskmax = itaskmin + ntask_per_rank_max - 1;
        }
        else
        {
            itaskmin = min_rank * ntask_per_rank_max + (r - min_rank) * ntask_per_rank_min;
            itaskmax = itaskmin + ntask_per_rank_min - 1;
        }

        for(int i = itaskmin, j = 0; i <= itaskmax; i++, j++)
            task_vector_partitioned[r][j] = i;
    }

    return task_vector_partitioned;
}

VVI util::general::partition_total_task_vector(VI task_vector, int world_size)
{
    int ntask = task_vector.size();
    VVI task_vector_partitioned(world_size);
    int ntask_per_rank_min = ntask / world_size;
    int ntask_per_rank_max = ntask % world_size == 0 ? ntask_per_rank_min : ntask_per_rank_min + 1;
    int min_rank = ntask % world_size;

    for(int r = 0; r < world_size; r++)
    {
        int ntask_my_rank = r < min_rank ? ntask_per_rank_max : ntask_per_rank_min;
        task_vector_partitioned[r].resize(ntask_my_rank);
        int itaskmin, itaskmax;

        if(r < min_rank)
        {
            itaskmin = r * ntask_per_rank_max;
            itaskmax = itaskmin + ntask_per_rank_max - 1;
        }
        else
        {
            itaskmin = min_rank * ntask_per_rank_max + (r - min_rank) * ntask_per_rank_min;
            itaskmax = itaskmin + ntask_per_rank_min - 1;
        }

        for(int i = itaskmin, j = 0; i <= itaskmax; i++, j++)
            task_vector_partitioned[r][j] = task_vector[i];

    }

    return task_vector_partitioned;
}
