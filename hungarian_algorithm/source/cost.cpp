// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include "../include/cost.hpp"

#include <algorithm> // std::nth_element
#include <functional> // std::function
#include <limits> // std::numeric_limits
#include <list> // std::list
#include <utility> // std::pair
#include <vector> // std::vector

using std::function;
using std::pair;
using std::vector;

typedef vector<int> VI;
typedef vector<double> VD;
typedef vector<VD> VVD;
typedef pair<int, double> PID;
typedef vector<PID> VPID;
typedef vector<VPID> VVPID;

// Maximum cost, symbolizing infinity
constexpr double cost_max = (double)std::numeric_limits<int>::max();

/* Compute cost objects */

VVD hungarian_algorithm::cost::create_costobject_adjmatrix(VVD &Xa, VVD &Xb, VD &p, function<double(VD&, VD&, VD&)> d)
{
    int N = Xa.size();

    // Define cost matrix object
    VVD cost_adjmatrix(N, VD(N, 0.0));

    // Compute cost matrix for each pair (i, j)
    // from the distance d
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            cost_adjmatrix[i][j] = d(Xa[i], Xb[j], p);
        }
    }

    return cost_adjmatrix;
}

VVD hungarian_algorithm::cost::create_costobject_adjmatrix(int N_edge, VVD &Xa, VVD &Xb, VD &p, function<double(VD&, VD&, VD&)> d)
{
    int N = Xa.size();

    // Define cost matrix object with standard value cost_max
    VVD cost_adjmatrix(N, VD(N, cost_max));

    // Compute all distances and set the N_edge smallest ones for each row
    VPID cost_node_row(N);
    // For each row i do:
    for(int i = 0; i < N; i++)
    {
        // For each column j do:
        // Compute the distance d_ij and store it together with the
        // corresponding index j into the j-th entry of cost_node_row
        for(int j = 0; j < N; j++)
        {
            PID pid;
            pid.first = j;
            pid.second = d(Xa[i], Xb[j], p);
            cost_node_row[j] = pid;
        }

        // Partial sort up to the n-th element of cost_node_row
        std::nth_element(
            cost_node_row.begin(), 
            cost_node_row.begin() + N_edge, 
            cost_node_row.end(), 
            [](PID pid_a, PID pid_b) {return pid_a.second < pid_b.second;}
        );

        // Assign N_edge smallest values in cost_adjmatrix
        for(int j = 0; j < N_edge; j++)
            cost_adjmatrix[i][cost_node_row[j].first] = cost_node_row[j].second;
    }

    // Do the same column-wise
    VPID cost_node_col(N);
    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < N; i++)
        {
            PID pid;
            pid.first = i;
            pid.second = d(Xa[i], Xb[j], p);
            cost_node_col[i] = pid;
        }

        std::nth_element(
            cost_node_col.begin(), 
            cost_node_col.begin() + N_edge, 
            cost_node_col.end(), 
            [](PID pid_a, PID pid_b) {return pid_a.second < pid_b.second;}
        );

        for(int i = 0; i < N_edge; i++)
            cost_adjmatrix[cost_node_col[i].first][j] = cost_node_col[i].second;
    }

    return cost_adjmatrix;
}

VVPID hungarian_algorithm::cost::create_costobject_adjlist_plain(int N_edge, VVD &Xa, VVD &Xb, VD &p, function<double(VD&, VD&, VD&)> d)
{
    int N = Xa.size();

    // Calculate
    VVD cost_adjmatrix = create_costobject_adjmatrix(N_edge, Xa, Xb, p, d);

    // Define cost matrix object with standard value cost_max
    VVPID cost_adjlist(N);

    int ctr;
    for(int i = 0; i < N; i++)
    {
        ctr = 0;
        for(int j = 0; j < N; j++)
        {
            ctr += cost_adjmatrix[i][j] < cost_max ? 1 : 0;
        }
        cost_adjlist[i].resize(ctr);

        for(int j = 0, k = 0; j < N; j++)
        {
            if(cost_adjmatrix[i][j] < cost_max)
            {
                cost_adjlist[i][k] = std::make_pair(j, cost_adjmatrix[i][j]);
                k++;
            }
        }
    }

    return cost_adjlist;
}

VVPID hungarian_algorithm::cost::create_costobject_adjlist_slim(int N_edge, VVD &Xa, VVD &Xb, VD &p, function<double(VD&, VD&, VD&)> d)
{
    int N = Xa.size();

    // Define cost matrix object with standard value cost_max
    VVPID cost_adjlist(N, VPID(2 * N_edge));

    // Find the N_edge nearest neighbours per row
    VVPID cost_adjlist_row(N, VPID(N_edge));
    VPID cost_node_row(N);
    VPID cost_node_row_part(N_edge);
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            cost_node_row[j] = std::make_pair(j, d(Xa[i], Xb[j], p));

        std::nth_element(cost_node_row.begin(), cost_node_row.begin() + N_edge, cost_node_row.end(), [](PID pid_a, PID pid_b) {return pid_a.second < pid_b.second;});
        std::copy(cost_node_row.begin(), cost_node_row.begin() + N_edge, cost_node_row_part.begin());
        std::sort(cost_node_row_part.begin(), cost_node_row_part.end(), [](PID pid_a, PID pid_b) {return pid_a.second < pid_b.second;});
        std::copy(cost_node_row_part.begin(), cost_node_row_part.end(), cost_adjlist_row[i].begin());
    }

    // Find the M nearest neighbours per column
    VVPID cost_adjlist_col(N, VPID(N_edge));
    VPID cost_node_col(N);
    VPID cost_node_col_part(N_edge);
    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < N; i++)
            cost_node_col[i] = std::make_pair(i, d(Xa[i], Xb[j], p));

        std::nth_element(cost_node_col.begin(), cost_node_col.begin() + N_edge, cost_node_col.end(), [](PID pid_a, PID pid_b) {return pid_a.second < pid_b.second;});
        std::copy(cost_node_col.begin(), cost_node_col.begin() + N_edge, cost_node_col_part.begin());
        std::sort(cost_node_col_part.begin(), cost_node_col_part.end(), [](PID pid_a, PID pid_b) {return pid_a.second < pid_b.second;});
        std::copy(cost_node_col_part.begin(), cost_node_col_part.end(), cost_adjlist_col[j].begin());
    }

    // Merge such that the number of neighbours is at least m row-&columnwise
    vector<std::list<PID>> cost_adjlist_tmp(N, std::list<PID>());
    for(int i = 0; i < N; i++)
        cost_adjlist_tmp[i].assign(cost_adjlist_row[i].begin(), cost_adjlist_row[i].end());

    bool found;
    for(int j = 0; j < N; j++)
    {
        for(VPID::iterator i_it = cost_adjlist_col[j].begin(); i_it != cost_adjlist_col[j].end(); i_it++)
        {
            int i = i_it->first;
            found = false;
            for(VPID::iterator j_it = cost_adjlist_row[i].begin(); j_it != cost_adjlist_row[i].end(); j_it++)
            {
                if(j == j_it->first)
                {
                    found = true;
                    break;
                }
            }
            if(!found)
                cost_adjlist_tmp[i].push_back(std::make_pair(j, i_it->second));
        }
    }

    // Sort the merged lists
    for(int i = 0; i < N; i++)
        cost_adjlist_tmp[i].sort([](PID pid_a, PID pid_b) {return pid_a.first < pid_b.first;});

    // Assign lists from cost_adjlist_tmp to vectors of cost_adjlist
    for(int i = 0; i < N; i++)
        cost_adjlist[i].assign(cost_adjlist_tmp[i].begin(), cost_adjlist_tmp[i].end());

    return cost_adjlist;
}

/* Compute total cost from given matching solution */

double hungarian_algorithm::cost::compute_total_cost_adjmatrix(VVD& cost_adjmatrix, VI q)
{
    int N = cost_adjmatrix.size();
    if(N != q.size())
        return std::numeric_limits<double>::max();

    double total_cost = 0.0;
    int count = 0;

    for(int i = 0; i < N; i++)
    {
        if(q[i] != -1)
        {
            count++;
            total_cost += cost_adjmatrix[i][q[i]];
        }
    }

    if(count != N)
        total_cost *= ((double) N / (double) count);

    return total_cost;
}

double hungarian_algorithm::cost::compute_total_cost_adjlist(VVPID& cost_adjlist, VI q)
{
    int N = cost_adjlist.size();
    if(N != q.size())
        return std::numeric_limits<double>::max();

    double total_cost = 0.0;
    int count = 0;
    int size_i;

    for(int i = 0; i < N; i++)
    {
        if(q[i] != -1)
        {
            size_i = cost_adjlist[i].size();
            bool found = false;
            for(int j = 0; j < size_i && !found; j++)
            {
                if(q[i] == cost_adjlist[i][j].first)
                {
                    total_cost += cost_adjlist[i][j].second;
                    count++;
                    found = true;
                }
            }
        }
    }

    if(count != N)
        total_cost *= ((double) N / (double) count);

    return total_cost;
}
