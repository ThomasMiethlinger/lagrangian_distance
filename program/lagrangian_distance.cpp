// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

// STL
#include <cmath> // sqrt
#include <iostream>  // std::cout, std::endl
#include <fstream> // std::ifstream
#include <string> // std::string
#include <functional> // std::function
#include <limits> // std::numeric_limits
#include <utility> // std::pair
#include <vector> // std::vector
#include <list> // std::list

// MPI
#include <mpi.h>

// Boost
#include <boost/algorithm/string.hpp>

// Internal classes
#include "../hungarian_algorithm/include/distance.hpp"
#include "../hungarian_algorithm/include/cost.hpp"
#include "../hungarian_algorithm/include/core.hpp"

#include "../util/include/io.hpp"
#include "../util/include/general.hpp"

using std::pair;
using std::size_t;
using std::string;
using std::vector;
using std::list;
using std::cout;
using std::endl;

typedef vector<int> VI;
typedef vector<VI> VVI;
typedef vector<double> VD;
typedef vector<VD> VVD;
typedef pair<int, double> PID;
typedef vector<PID> VPID;
typedef vector<VPID> VVPID;

int main(int argc, char * argv[])
{
    // Initialize MPI environment and query for rank and size
    MPI_Init(&argc, &argv);
    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Variables for the processing of the input dataset
    // They are set by the programm arguments argv
    int N; // Number of particles
    int N_vertex; // Reduced number of particles
    double ratio_N;
    int N_edge; // Number of particles to consider for matching
    VD p(6);
    std::function<double(VD&, VD&, VD&)> d;
    int tmin;
    int tmax;
    int tstep;
    int nsteps;
    int ntask;
    vector<pair<int, int>> total_task_vector; // Define a vector total_task_vector which transforms from i to {t_1, t_2} for all tasks i

    // IO
    string input_folder_path;
    string input_file_name_part;
    string input_file_path_name;

    bool output_matrix;
    string output_file_name_matrix;
    string output_folder_path_matrix;
    string output_file_path_name_matrix;
    bool output_list;
    string output_file_name_list;
    string output_folder_path_list;
    string output_file_path_name_list;
    bool output_matching;
    string output_file_name_matching;
    string output_folder_path_matching;
    string output_file_path_name_matching;

    // Read config file
    vector<string> argsv;
    argsv.assign(argv, argv + argc);
    string input_config_file = argsv[1];

    std::ifstream filestream;
    filestream.open(input_config_file);

    if(filestream.is_open())
    {
        string line;
        vector<string> parts;
        for (int i = 0; std::getline(filestream, line); i++)
        {
            boost::split(parts, line, boost::is_any_of(" "));
            if(parts[0] == "numbers") // numbers 50000 100 100
            {
                N = std::stoi(parts[1]);
                N_vertex = std::stoi(parts[2]);
                ratio_N = (double)N / (double)N_vertex;
                N_edge = std::stoi(parts[3]);
            }
            else if(parts[0] == "distance") // distance l2_distance_weight 1,1,1,0,0,0
            {
                d = hungarian_algorithm::distance::get_distance(parts[1]);
                vector<string> p_vec;
                boost::split(p_vec, parts[2], boost::is_any_of(","));
                std::transform(p_vec.begin(), p_vec.end(), p.begin(), [](string s) { return std::stod(s); });
            }
            else if(parts[0] == "time") // time matrix 10000 20000 2000 OR time list time_pairs.txt
            {
                if(parts[1] == "matrix") // time matrix 10000 20000 2000
                {
                    tmin = std::stoi(parts[2]);
                    tmax = std::stoi(parts[3]);
                    tstep = std::stoi(parts[4]);
                    nsteps = (tmax - tmin) / tstep + 1;
                    ntask = (nsteps * nsteps - nsteps) / 2;
                    total_task_vector.resize(ntask);
                    for(int t1 = tmin, i = 0; t1 <= (tmax - tstep); t1 += tstep)
                    {
                        for(int t2 = t1 + tstep; t2 <= tmax; t2 += tstep, i++)
                        {
                            total_task_vector[i] = std::make_pair(t1, t2);
                        }
                    }
                }
                else if(parts[1] == "list") // time list time_pairs.txt
                {
                    list<pair<int, int>> time_pair_list = util::io::read_total_task_time_pairs(parts[2]);
                    total_task_vector.assign(time_pair_list.begin(), time_pair_list.end());
                }
            }
            else if(parts[0] == "input") // input dump% ../data/test_case/
            {
                input_file_name_part = parts[1];
                input_folder_path = parts[2];
            }
            else if(parts[0] == "matrix") // matrix 1 lagrangian_distance_matrix.txt ../results/distance_matrix/
            {
                output_matrix = std::stoi(parts[1].c_str());
                output_file_name_matrix = parts[2];
                output_folder_path_matrix = parts[3];
                output_file_path_name_matrix = output_folder_path_matrix + output_file_name_matrix;
            }
            else if(parts[0] == "list") // list 1 lagrangian_distance_%.txt ../results/distance_list/
            {
                output_list = std::stoi(parts[1].c_str());
                output_file_name_list = parts[2].replace(parts[2].find("%"), 1, std::to_string(world_rank));
                output_folder_path_list = parts[3];
                output_file_path_name_list = output_folder_path_list + output_file_name_list;
            }
            else if(parts[0] == "matching") // matching 1 matching_%.txt ../results/distance_matching/
            {
                output_matching = std::stoi(parts[1].c_str());
                output_file_name_matching = parts[2];
                output_folder_path_matching = parts[3];
            }
        }
        filestream.close();
    }
    else
    {
        cout << "Error: Could not read file " << input_config_file << "!" << endl;
        
        MPI_Finalize();
        return 1;
    }

    // Let rank 0 create the output folders, if they do not exist
    if(world_rank == 0)
    {
        util::io::create_work_directory(output_folder_path_matrix.c_str());
        util::io::create_work_directory(output_folder_path_list.c_str());
        util::io::create_work_directory(output_folder_path_matching.c_str());
    }


    // Vector for tasking
    VI task_vector_my_rank;
    int ntask_my_rank;

    // Rank 0 needs to read the intermediate results, i.e., look up existing results whats there
    if(world_rank == 0)
    {
        list<pair<pair<int, int>, double>> total_ready_results;
        if(!util::io::get_total_ready_results(total_ready_results, output_folder_path_list))
        {
            MPI_Finalize();
            return 1;
        }

        list<pair<int, int>> total_ready_time_pairs;
        for(auto it = total_ready_results.begin(); it != total_ready_results.end(); it++)
            total_ready_time_pairs.push_back((*it).first);

        VI total_task_indices(ntask);
        int ctr_ready_results = util::general::mark_ready_task(total_task_vector, total_ready_time_pairs, total_task_indices);

        VI total_remaining_task_indices(ntask - ctr_ready_results);

        for(int i = 0, j = 0; i < ntask; i++)
        {
            if(total_task_indices[i] != -1)
            {
                total_remaining_task_indices[j] = total_task_indices[i];
                j++;
            }
        }
        VVI task_index_vector_partitioned = util::general::partition_total_task_vector(total_remaining_task_indices, world_size);

        task_vector_my_rank = task_index_vector_partitioned[0];
        ntask_my_rank = task_vector_my_rank.size();
        for(int r = 1; r < world_size; r++)
        {
            int ntask_r = task_index_vector_partitioned[r].size();
            MPI_Ssend(&ntask_r, 1, MPI_INT, r, 0, MPI_COMM_WORLD);
            MPI_Ssend(task_index_vector_partitioned[r].data(), task_index_vector_partitioned[r].size(), MPI_INT, r, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&ntask_my_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        task_vector_my_rank.resize(ntask_my_rank);
        MPI_Recv(task_vector_my_rank.data(), ntask_my_rank, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // For each rank and its task, save the result of the computation into lagrangian_distance_results
    VD lagrangian_distance_results(ntask_my_rank, 0);

    // Initialize particle state vectors Xa and Xb
    VVD Xa(N, VD(6));
    VVD Xb(N, VD(6));
    VVD Xa_vertex(N_vertex, VD(6));
    VVD Xb_vertex(N_vertex, VD(6));

    // Variables for the indices of each task
    int t1, t2;
    int t1_prev, t2_prev;
    t1_prev = t2_prev = std::numeric_limits<int>::max();
    pair<int, int> task_time_pair;
    pair<pair<int, int>, double> result_pair;

    for(int i = 0; i < ntask_my_rank; i++)
    {
        // New task: computing lagrangian distance of (t_1, t_2) = task_time_pair
        task_time_pair = total_task_vector[task_vector_my_rank[i]];
        
        t1 = task_time_pair.first;
        string s1 = input_file_name_part;
        // In order to avoid two-times reading from the same file
        if(t1 != t1_prev)
        {
            input_file_path_name = input_folder_path + s1.replace(s1.find("%"), 1, std::to_string(t1));
            // Read input dump file
            if(!util::io::read_state(input_file_path_name, Xa))
            {
                MPI_Finalize();
                return 1;
            }
            for(int i = 0; i < N_vertex; i++)
                Xa_vertex[i] = Xa[(int)(i * ratio_N)];
        }

        // Same for state b
        t2 = task_time_pair.second;
        string s2 = input_file_name_part;
        if(t2 != t2_prev)
        {
            input_file_path_name = input_folder_path + s2.replace(s2.find("%"), 1, std::to_string(t2));
            if(!util::io::read_state(input_file_path_name, Xb))
            {
                MPI_Finalize();
                return 1;
            }
            for(int i = 0; i < N_vertex; i++)
	            Xb_vertex[i] = Xb[(int)(i * ratio_N)];
        }

        VI q;
        if(N_vertex == N_edge)
        {
            // Compute the distance adjacency lists
            VVD cost_adjmatrix = hungarian_algorithm::cost::create_costobject_adjmatrix(Xa_vertex, Xb_vertex, p, d);

            // Compute the hungarian algorithm matching solution for adjacency lists
            q = hungarian_algorithm::core::min_cost_matching_adjmatrix(cost_adjmatrix);

            // Compute the total cost from the matching solution
            lagrangian_distance_results[i] = hungarian_algorithm::cost::compute_total_cost_adjmatrix(cost_adjmatrix, q);


            if(output_matching)
            {
                string sm = output_file_name_matching;
                output_file_path_name_matching = output_folder_path_matching + sm.replace(sm.find("%"), 1, std::to_string(t1) + "_" + std::to_string(t2));
                util::io::write_matching_result(output_file_path_name_matching, Xa_vertex, Xb_vertex, cost_adjmatrix, lagrangian_distance_results[i], q);
            }
        }
        else
        {
            // Compute the distance adjacency lists
            VVPID cost_adjlist = hungarian_algorithm::cost::create_costobject_adjlist_plain(N_edge, Xa_vertex, Xb_vertex, p, d);

            // Compute the hungarian algorithm matching solution for adjacency lists
            q = hungarian_algorithm::core::min_cost_matching_adjlist(cost_adjlist);

            // Compute the total cost from the matching solution
            lagrangian_distance_results[i] = hungarian_algorithm::cost::compute_total_cost_adjlist(cost_adjlist, q);
        }

        // Write partial results
        result_pair.first = task_time_pair;
        result_pair.second = lagrangian_distance_results[i];

        if(output_list)
            util::io::write_list_result(output_file_path_name_list, result_pair);

        // Update indices
        t1_prev = t1;
        t2_prev = t2;
    }

    // If a full distance matrix should be written as output
    if(output_matrix)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        // Rank 0 needs to compose the lagrangian_distance_matrix
        if(world_rank == 0)
        {
            list<pair<pair<int, int>, double>> lagrangian_distance_results;
            if(!util::io::get_total_ready_results(lagrangian_distance_results, output_folder_path_list))
            {
                MPI_Finalize();
                return 1;
            }

            // Set the final lagrangian_distance_matrix from received results
            VVD lagrangian_distance_matrix(nsteps, VD(nsteps, 0));
            pair<pair<int, int>, double> task;
            pair<int, int> task_time_pair;
            int i1, i2;
            double task_res;
            for(list<pair<pair<int, int>, double>>::iterator it = lagrangian_distance_results.begin(); it != lagrangian_distance_results.end(); it++)
            {
                task = *it;
                task_time_pair = task.first;
                i1 = (task_time_pair.first - tmin) / tstep;
                i2 = (task_time_pair.second - tmin) / tstep;
                task_res = task.second;
                lagrangian_distance_matrix[i1][i2] = lagrangian_distance_matrix[i2][i1] = task_res;
            }
            // Write the computed distance matrix based on discretes into this filepathname
            util::io::write_matrix_result(output_file_path_name_matrix, lagrangian_distance_matrix);
        }
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}
