// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include <utility> // std::pair
#include <vector> // std::vector

#pragma once

namespace hungarian_algorithm
{
    namespace core
    {
        std::vector<int> min_cost_matching_adjmatrix(
            std::vector<std::vector<double>> &cost_adjmatrix
        );

        std::vector<int> min_cost_matching_adjlist(
            std::vector<std::vector<std::pair<int, double>>> &cost_adjlist
        );
    }
};
