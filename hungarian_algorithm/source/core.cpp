// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include "../include/core.hpp"

#include <cmath> // std::fabs
#include <functional> // std::function
#include <list> // std::list
#include <limits> // std::numeric_limits
#include <set> // std::set
#include <utility> // std::pair
#include <vector> // std::vector

using std::function;
using std::list;
using std::pair;
using std::set;
using std::vector;

typedef vector<int> VI;
typedef vector<double> VD;
typedef vector<VD> VVD;
typedef list<int> LI;
typedef set<int> SI;
typedef pair<int, double> PID;
typedef vector<PID> VPID;
typedef vector<VPID> VVPID;

// Precision which is used to check for zero cost
constexpr double precision = 0.000000001;

// Maximum cost, symbolizing infinity
constexpr double cost_max = (double)std::numeric_limits<int>::max();

// Compute the hungarian algorithm minimum cost matching solution
// for an adjacency cost matrix
VI hungarian_algorithm::core::min_cost_matching_adjmatrix(VVD &cost_adjmatrix)
{
    int N = cost_adjmatrix.size();

    // construct dual feasible solution
    VD u(N, cost_max);
    VD v(N, cost_max);
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            // Minimize cost_adjmatrix for each row i w.r.t. all columns j
            u[i] = std::min(u[i], cost_adjmatrix[i][j]);
        }
    }
    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < N; i++)
        {
            // Minimize cost_adjmatrix - u for each column j w.r.t. all rows i
            v[j] = std::min(v[j], cost_adjmatrix[i][j] - u[i]);
        }
    }

    // construct primal solution satisfying complementary slackness
    VI q = VI(N, -1);
    VI q_inv = VI(N, -1);
    int mated = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            // if q_inv[j] has a matching we can continue
            if(q_inv[j] != -1)
                continue;
            // if cost_adjmatrix[i][j] - u[i] - v[j] == 0
            // we match (i, j) and increase the counter mated
            if(fabs(cost_adjmatrix[i][j] - u[i] - v[j]) < precision)
            {
                q[i] = j;
                q_inv[j] = i;
                mated++;
                break;
            }
        }
    }

    VD dist(N);
    VI dad(N);
    VI seen(N);

    // repeat until primal solution is feasible
    while(mated < N)
    {
        // find an unmatched left node
        // store the corresponding index of this single node in s
        int s = 0;
        while(q[s] != -1)
            s++;

        // initialize Dijkstra
        std::fill(dad.begin(), dad.end(), -1);
        std::fill(seen.begin(), seen.end(), 0);
        for(int k = 0; k < N; k++)
        {
            dist[k] = cost_adjmatrix[s][k] - u[s] - v[k];
        }

        int j = 0;

        while(true)
        {
            // find closest
            j = -1;

            for(int k = 0; k < N; k++)
            {
                if(!seen[k])
                {
                    if(j == -1 || dist[k] < dist[j])
                        j = k;
                }
            }
            seen[j] = 1;

            // termination condition
            if(q_inv[j] == -1)
                break;

            // relax neighbours
            const int i = q_inv[j];
            for(int k = 0; k < N; k++)
            {
                if(!seen[k])
                {
                    const double new_dist = dist[j] + cost_adjmatrix[i][k] - u[i] - v[k];
                    if(new_dist < dist[k])
                    {
                        dist[k] = new_dist;
                        dad[k] = j;
                    }
                }
            }
        }

        // update dual variables
        for(int k = 0; k < N; k++)
        {
            if(k != j && seen[k])
            {
                const double delta = dist[k] - dist[j];
                v[k] += delta;
                u[q_inv[k]] -= delta;
            }
        }
        u[s] += dist[j];

        // augment along path
        while(dad[j] >= 0)
        {
            const int d = dad[j];
            q_inv[j] = q_inv[d];
            q[q_inv[j]] = j;
            j = d;
        }
        q_inv[j] = s;
        q[s] = j;

        mated++;
    }

    return q;
}

VI min_cost_matching_adjlist_backward(VVPID &cost_adjlist)
{
    int N = cost_adjlist.size();

    VD u(N, cost_max);
    VD v(N, cost_max);
    for(int i = 0; i < N; i++)
    {
        const int size = cost_adjlist[i].size();
        for(int j = 0; j < size; j++)
        {
            u[i] = std::min(u[i], cost_adjlist[i][j].second);
        }
    }
    for(int i = 0; i < N; i++)
    {
        const int size = cost_adjlist[i].size();
        for(int j = 0, y; j < size; j++)
        {
            y = cost_adjlist[i][j].first;
            v[y] = std::min(v[y], cost_adjlist[i][j].second - u[i]);
        }
    }

    VI q = VI(N, -1);
    VI q_inv = VI(N, -1);
    int mated = 0;
    for(int i = 0; i < N; i++)
    {
        const int size = cost_adjlist[i].size();
        for(int j = 0, y; j < size; j++)
        {
            y = cost_adjlist[i][j].first;
            if(q_inv[y] != -1)
                continue;

            if(fabs(cost_adjlist[i][j].second - u[i] - v[y]) < precision)
            {
                q[i] = y;
                q_inv[y] = i;
                mated++;
                break;
            }
        }
    }

    VD dist(N);
    VI dad(N);
    VI seen(N);
    LI left_vert;

    while(mated < N)
    {
        left_vert.clear();

        int s = 0;
        while(q[s] != -1)
            s++;

        const int size_s = cost_adjlist[s].size();
        std::fill(dad.begin(), dad.end(), -1);
        std::fill(seen.begin(), seen.end(), 0);
        std::fill(dist.begin(), dist.end(), cost_max);
        for(int k = 0, h; k < size_s; k++)
        {
            h = cost_adjlist[s][k].first;
            dist[h] = cost_adjlist[s][k].second - u[s] - v[h];
        }

        int j = 0;
        int i = -1;
        while(true)
        {
            j = -1;

            if(i == -1)
                i = s;                
            left_vert.push_front(i);

            for(LI::iterator it = left_vert.begin(); it != left_vert.end() && j == -1; it++)
            {
                int i_past = *it;
                int size_i_past = cost_adjlist[i_past].size();

                for(int k = 0, h; k < size_i_past; k++)
                {
                    h = cost_adjlist[i_past][k].first;

                    if(!seen[h])
                    {
                        if(j == -1 || dist[h] < dist[j])
                            j = h;
                    }
                }
            }

            if(j == -1)
                return q;

            seen[j] = 1;

            if(q_inv[j] == -1)
                break;
            i = q_inv[j];
            int size_i = cost_adjlist[i].size();
            for(int k = 0, h; k < size_i; k++)
            {
                h = cost_adjlist[i][k].first;
                if(!seen[h])
                {
                    const double new_dist = dist[j] + cost_adjlist[i][k].second - u[i] - v[h];
                    if(new_dist < dist[h])
                    {
                        dist[h] = new_dist;
                        dad[h] = j;
                    }
                }
            }
        }

        for(int h = 0; h < N; h++)
        {
            if(h != j && seen[h])
            {
                const double delta = dist[h] - dist[j];
                v[h] += delta;
                u[q_inv[h]] -= delta;
            }
        }
        u[s] += dist[j];

        while(dad[j] >= 0)
        {
            const int d = dad[j];
            q_inv[j] = q_inv[d];
            q[q_inv[j]] = j;
            j = d;
        }
        q_inv[j] = s;
        q[s] = j;

        mated++;
    }

    return q;
}


VI min_cost_matching_adjlist_minset(int N, VVPID &cost_adjlist)
{
    function<bool(PID, PID)> comparator_PID = [](PID x1, PID x2){ return x1.second != x2.second ? (x1.second < x2.second) : (x1.first < x2.first); };

    VD u(N, cost_max);
    VD v(N, cost_max);
    for(int i = 0; i < N; i++)
    {
        const int size = cost_adjlist[i].size();
        for(int j = 0; j < size; j++)
        {
            u[i] = std::min(u[i], cost_adjlist[i][j].second);
        }
    }
    for(int i = 0; i < N; i++)
    {
        const int size = cost_adjlist[i].size();
        for(int j = 0, y; j < size; j++)
        {
            y = cost_adjlist[i][j].first;
            v[y] = std::min(v[y], cost_adjlist[i][j].second - u[i]);
        }
    }

    VI q = VI(N, -1);
    VI q_inv = VI(N, -1);
    int mated = 0;
    for(int i = 0; i < N; i++)
    {
        const int size = cost_adjlist[i].size();
        for(int j = 0, y; j < size; j++)
        {
            y = cost_adjlist[i][j].first;
            if(q_inv[y] != -1)
                continue;

            if(fabs(cost_adjlist[i][j].second - u[i] - v[y]) < precision)
            {
                q[i] = y;
                q_inv[y] = i;
                mated++;
                break;
            }
        }
    }

    VD dist(N);
    VI dad(N);
    VI seen(N);

    set<PID, function<bool(PID, PID)>> dist_cache_ordered(comparator_PID);
    pair<bool, set<PID, function<bool(PID, PID)>>::iterator> inactive_pair = std::make_pair(false, dist_cache_ordered.end());
    vector<pair<bool, set<PID, function<bool(PID, PID)>>::iterator>> dist_active(N, inactive_pair);
    while(mated < N)
    {
        dist_cache_ordered.clear();
        std::fill(dist_active.begin(), dist_active.end(), inactive_pair);

        int s = 0;
        while(q[s] != -1)
            s++;

        const int size_s = cost_adjlist[s].size();
        std::fill(dad.begin(), dad.end(), -1);
        std::fill(seen.begin(), seen.end(), 0);
        std::fill(dist.begin(), dist.end(), cost_max);
        for(int k = 0, h; k < size_s; k++)
        {
            h = cost_adjlist[s][k].first;
            dist[h] = cost_adjlist[s][k].second - u[s] - v[h];
            if(dist_active[h].first)
            {
                dist_cache_ordered.erase(dist_active[h].second);
                dist_active[h].second = dist_cache_ordered.insert(std::make_pair(h, dist[h])).first;
            }
            else
            {
                dist_active[h].first = true;
                dist_active[h].second = dist_cache_ordered.insert(std::make_pair(h, dist[h])).first;
            }
        }

        int j;
        int i = -1;
        while(true)
        {
            j = -1;

            if(i == -1)
                i = s;          

            set<PID>::iterator it_dist_min = dist_cache_ordered.begin();
            if(it_dist_min != dist_cache_ordered.end())
            {
                j = it_dist_min->first;
                dist_active[j].first = false;
                dist_active[j].second = dist_cache_ordered.end();
                dist_cache_ordered.erase(it_dist_min);
            }

            if(j == -1)
                return q;

            seen[j] = 1;

            if(q_inv[j] == -1)
                break;

            i = q_inv[j];
            int size = cost_adjlist[i].size();
            for(int k = 0, h; k < size; k++)
            {
                h = cost_adjlist[i][k].first;
                if(!seen[h])
                {
                    const double new_dist = dist[j] + cost_adjlist[i][k].second - u[i] - v[h];
                    if(new_dist < dist[h])
                    {
                        dist[h] = new_dist;
                        dad[h] = j;

                        if(dist_active[h].first)
                        {
                            dist_cache_ordered.erase(dist_active[h].second);
                            dist_active[h].second = dist_cache_ordered.insert(std::make_pair(h, dist[h])).first;
                        }
                        else
                        {
                            dist_active[h].first = true;
                            dist_active[h].second = dist_cache_ordered.insert(std::make_pair(h, dist[h])).first;
                        }
                    }
                }
            }
        }

        for(int h = 0; h < N; h++)
        {
            if(h != j && seen[h])
            {
                const double delta = dist[h] - dist[j];
                v[h] += delta;
                u[q_inv[h]] -= delta;
            }
        }
        u[s] += dist[j];

        while(dad[j] >= 0)
        {
            const int d = dad[j];
            q_inv[j] = q_inv[d];
            q[q_inv[j]] = j;
            j = d;
        }
        q_inv[j] = s;
        q[s] = j;

        mated++;
    }

    return q;
}

VI min_cost_matching_adjlist_minset_limiter(int N, VVPID &cost_adjlist, int size_limit_cache)
{
    function<bool(PID, PID)> comparator_PID = [](PID x1, PID x2){ return x1.second != x2.second ? (x1.second < x2.second) : (x1.first < x2.first); };

    VD u(N, cost_max);
    VD v(N, cost_max);
    for(int i = 0; i < N; i++)
    {
        const int size = cost_adjlist[i].size();
        for(int j = 0; j < size; j++)
        {
            u[i] = std::min(u[i], cost_adjlist[i][j].second);
        }
    }
    for(int i = 0; i < N; i++)
    {
        const int size = cost_adjlist[i].size();
        for(int j = 0, y; j < size; j++)
        {
            y = cost_adjlist[i][j].first;
            v[y] = std::min(v[y], cost_adjlist[i][j].second - u[i]);
        }
    }

    VI q = VI(N, -1);
    VI q_inv = VI(N, -1);
    int mated = 0;
    for(int i = 0; i < N; i++)
    {
        const int size = cost_adjlist[i].size();
        for(int j = 0, y; j < size; j++)
        {
            y = cost_adjlist[i][j].first;
            if(q_inv[y] != -1)
                continue;

            if(fabs(cost_adjlist[i][j].second - u[i] - v[y]) < precision)
            {
                q[i] = y;
                q_inv[y] = i;
                mated++;
                break;
            }
        }
    }

    VD dist(N);
    VI dad(N);
    VI seen(N);

    set<PID, function<bool(PID, PID)>> dist_cache_ordered(comparator_PID);
    pair<bool, set<PID, function<bool(PID, PID)>>::iterator> inactive_pair = std::make_pair(false, dist_cache_ordered.end());
    vector<pair<bool, set<PID, function<bool(PID, PID)>>::iterator>> dist_active(N, inactive_pair);
    set<PID>::iterator it_dist_max = dist_cache_ordered.end();
    while(mated < N)
    {
        dist_cache_ordered.clear();
        std::fill(dist_active.begin(), dist_active.end(), inactive_pair);

        int s = 0;
        while(q[s] != -1)
            s++;

        const int size_s = cost_adjlist[s].size();
        std::fill(dad.begin(), dad.end(), -1);
        std::fill(seen.begin(), seen.end(), 0);
        std::fill(dist.begin(), dist.end(), cost_max);
        for(int k = 0, h; k < size_s; k++)
        {
            h = cost_adjlist[s][k].first;
            dist[h] = cost_adjlist[s][k].second - u[s] - v[h];
            if(dist_active[h].first)
            {
                dist_cache_ordered.erase(dist_active[h].second);
                dist_active[h].second = dist_cache_ordered.insert(std::make_pair(h, dist[h])).first;
            }
            else
            {
                dist_active[h].first = true;
                dist_active[h].second = dist_cache_ordered.insert(std::make_pair(h, dist[h])).first;
            }
        }

        int j;
        int i = -1;
        while(true)
        {
            j = -1;

            if(i == -1)
                i = s;          

            set<PID>::iterator it_dist_min = dist_cache_ordered.begin();
            if(it_dist_min != dist_cache_ordered.end())
            {
                j = it_dist_min->first;
                dist_active[j].first = false;
                dist_active[j].second = dist_cache_ordered.end();
                dist_cache_ordered.erase(it_dist_min);
            }

            if(j == -1)
                return q;

            seen[j] = 1;

            if(q_inv[j] == -1)
                break;

            i = q_inv[j];
            int size = cost_adjlist[i].size();
            for(int k = 0, h; k < size; k++)
            {
                h = cost_adjlist[i][k].first;
                if(!seen[h])
                {
                    const double new_dist = dist[j] + cost_adjlist[i][k].second - u[i] - v[h];
                    if(new_dist < dist[h])
                    {
                        dist[h] = new_dist;
                        dad[h] = j;

                        if(dist_active[h].first)
                        {
                            dist_cache_ordered.erase(dist_active[h].second);
                            dist_active[h].second = dist_cache_ordered.insert(std::make_pair(h, dist[h])).first;
                        }
                        else
                        {
                            dist_active[h].first = true;
                            dist_active[h].second = dist_cache_ordered.insert(std::make_pair(h, dist[h])).first;
                        }

                        while((int)dist_cache_ordered.size() > size_limit_cache)
                        {
                            it_dist_max = dist_cache_ordered.end();
                            it_dist_max--;

                            int l = it_dist_max->first;
                            dist_active[l].first = false;
                            dist_active[l].second = dist_cache_ordered.end();
                            dist_cache_ordered.erase(it_dist_max);
                        }      
                    }
                }
            }
        }

        for(int h = 0; h < N; h++)
        {
            if(h != j && seen[h])
            {
                const double delta = dist[h] - dist[j];
                v[h] += delta;
                u[q_inv[h]] -= delta;
            }
        }
        u[s] += dist[j];

        while(dad[j] >= 0)
        {
            const int d = dad[j];
            q_inv[j] = q_inv[d];
            q[q_inv[j]] = j;
            j = d;
        }
        q_inv[j] = s;
        q[s] = j;

        mated++;
    }

    return q;
}

VI hungarian_algorithm::core::min_cost_matching_adjlist(VVPID &cost_adjlist)
{
    return min_cost_matching_adjlist_backward(cost_adjlist);
}
