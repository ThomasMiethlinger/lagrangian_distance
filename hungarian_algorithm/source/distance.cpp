// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include "../include/distance.hpp"

#include <algorithm> // std::transform
#include <functional> // std::function, std::plus, std::minus, std::multiplies
#include <cmath> // std::sqrt, std::fabs
#include <limits> // std::numeric_limits
#include <numeric> // std::accumulate
#include <string> // std::string

using std::vector;

typedef vector<double> VD;

std::function<double(VD&, VD&, VD&)> hungarian_algorithm::distance::get_distance(std::string const& s)
{
    if (s == "l1_distance_weight") return hungarian_algorithm::distance::l1_distance_weight;
    else if (s == "l2_distance_weight") return hungarian_algorithm::distance::l2_distance_weight;
    return nullptr;
}

/* Standard Distance Functions */

/*double hungarian_algorithm::distance::l1_distance(VD &xa, VD &xb)
{
    int na = xa.size();    
    int nb = xb.size();
    if(na != nb)
        return std::numeric_limits<double>::max();

    VD x(na);
    std::transform(xa.begin(), xa.end(), xb.begin(), x.begin(), std::minus<double>());
    std::transform(x.begin(), x.end(), x.begin(), fabs);
    return std::accumulate(x.begin(), x.end(), 0.0, std::plus<double>());
}

double hungarian_algorithm::distance::l2_distance(VD &xa, VD &xb)
{
    int na = xa.size();    
    int nb = xb.size();
    if(na != nb)
        return std::numeric_limits<double>::max();

    VD x(na);
    std::transform(xa.begin(), xa.end(), xb.begin(), x.begin(), std::minus<double>());
    std::transform(x.begin(), x.end(), x.begin(), x.begin(), std::multiplies<double>());
    return sqrt(std::accumulate(x.begin(), x.end(), 0.0, std::plus<double>()));
}*/

/*double hungarian_algorithm::distance::squared_l2_distance(VD &xa, VD &xb, VD &p)
{
    int na = xa.size();    
    int nb = xb.size();
    if(na != nb)
        return std::numeric_limits<double>::max();

    VD x(na);
    std::transform(xa.begin(), xa.end(), xb.begin(), x.begin(), std::minus<double>());
    std::transform(x.begin(), x.end(), x.begin(), square);
    std::transform(x.begin(), x.end(), p.begin(), x.begin(), std::multiplies<double>());
    return std::accumulate(x.begin(), x.end(), 0.0, std::plus<double>());
}*/

/* Standard Distance Functions with additional weights for each component */

double hungarian_algorithm::distance::l1_distance_weight(VD &xa, VD &xb, VD &p)
{
    int na = xa.size();    
    int nb = xb.size();
    if(na != nb)
        return std::numeric_limits<double>::max();

    VD x(na);
    std::transform(xa.begin(), xa.end(), xb.begin(), x.begin(), std::minus<double>());
    std::transform(x.begin(), x.end(), x.begin(), fabs);
    std::transform(x.begin(), x.end(), p.begin(), x.begin(), std::multiplies<double>());
    return std::accumulate(x.begin(), x.end(), 0.0, std::plus<double>());
}

double hungarian_algorithm::distance::l2_distance_weight(VD &xa, VD &xb, VD &p)
{
    int na = xa.size();    
    int nb = xb.size();
    if(na != nb)
        return std::numeric_limits<double>::max();

    VD x(na);
    std::transform(xa.begin(), xa.end(), xb.begin(), x.begin(), std::minus<double>());
    std::transform(x.begin(), x.end(), x.begin(), x.begin(), std::multiplies<double>());
    std::transform(x.begin(), x.end(), p.begin(), x.begin(), std::multiplies<double>());
    return sqrt(std::accumulate(x.begin(), x.end(), 0.0, std::plus<double>()));
}

double hungarian_algorithm::distance::l2_distance_squared_weight(VD &xa, VD &xb, VD &p)
{
    int na = xa.size();    
    int nb = xb.size();
    if(na != nb)
        return std::numeric_limits<double>::max();

    VD x(na);
    std::transform(xa.begin(), xa.end(), xb.begin(), x.begin(), std::minus<double>());
    std::transform(x.begin(), x.end(), x.begin(), x.begin(), std::multiplies<double>());
    std::transform(x.begin(), x.end(), p.begin(), x.begin(), std::multiplies<double>());
    return std::accumulate(x.begin(), x.end(), 0.0, std::plus<double>());
}
