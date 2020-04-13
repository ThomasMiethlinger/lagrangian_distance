// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include <functional> // std::function
#include <vector> // std::vector

#pragma once

namespace hungarian_algorithm
{
    namespace distance
    {
        std::function<double(std::vector<double>&, std::vector<double>&, std::vector<double>&)> get_distance(
            std::string const& s
        );

        /* Standard Distance Functions */
        /*double l1_distance(
            std::vector<double> &xa, 
            std::vector<double> &xb
        );

        double l2_distance(
            std::vector<double> &xa, 
            std::vector<double> &xb
        );*/

        /*double l2_distance_squared(
            std::vector<double> &xa, 
            std::vector<double> &xb
        );*/

        /*double lmax_distance(
            std::vector<double> &xa, 
            std::vector<double> &xb
        );*/

        /* Standard Distance Functions with additional weights for each component */

        double l1_distance_weight(
            std::vector<double> &xa, 
            std::vector<double> &xb, 
            std::vector<double> &p
        );

        double l2_distance_weight(
            std::vector<double> &xa, 
            std::vector<double> &xb, 
            std::vector<double> &p
        );

        double l2_distance_squared_weight(
            std::vector<double> &xa, 
            std::vector<double> &xb, 
            std::vector<double> &p
        );

        /*double squared_l2_distance_weight(
            std::vector<double> &xa, 
            std::vector<double> &xb, 
            std::vector<double> &p
        );*/

        /*double lmax_distance_weight(
            std::vector<double> &xa, 
            std::vector<double> &xb, 
            std::vector<double> &p
        );*/
    }
};

