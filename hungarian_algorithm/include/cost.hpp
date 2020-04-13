// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include <list> // std::list
#include <functional> // std::function
#include <utility> // std::pair
#include <vector> // std::vector

#pragma once

namespace hungarian_algorithm
{
    namespace cost
    {
        /* Compute cost objects */

        std::vector<std::vector<double>> create_costobject_adjmatrix(
            std::vector<std::vector<double>> &Xa, 
            std::vector<std::vector<double>> &Xb, 
            std::vector<double> &p, 
            std::function<double(std::vector<double>&, std::vector<double>&, std::vector<double>&)> d
        );

        std::vector<std::vector<double>> create_costobject_adjmatrix(
            int N_edge, 
            std::vector<std::vector<double>> &Xa, 
            std::vector<std::vector<double>> &Xb, 
            std::vector<double> &p, 
            std::function<double(std::vector<double>&, std::vector<double>&, std::vector<double>&)> d
        );

        std::vector<std::vector<std::pair<int, double>>> create_costobject_adjlist_plain(
            int N_edge, 
            std::vector<std::vector<double>> &Xa, 
            std::vector<std::vector<double>> &Xb, 
            std::vector<double> &p, 
            std::function<double(std::vector<double>&, std::vector<double>&, std::vector<double>&)> d
        );

        std::vector<std::vector<std::pair<int, double>>> create_costobject_adjlist_slim(
            int N_edge, 
            std::vector<std::vector<double>> &Xa, 
            std::vector<std::vector<double>> &Xb, 
            std::vector<double> &p, 
            std::function<double(std::vector<double>&, 
            std::vector<double>&, std::vector<double>&)> d
        );

        /* Compute total cost from given matching solution */

        double compute_total_cost_adjmatrix(
            std::vector<std::vector<double>>& cost_adjmatrix, 
            std::vector<int> q
        );

        double compute_total_cost_adjlist(
            std::vector<std::vector<std::pair<int, double>>>& cost_adjlist, 
            std::vector<int> q
        );
    }
};
