// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include <list> // std::list
#include <string> // std::string
#include <utility> // std::pair
#include <vector> // std::vector

#pragma once

namespace util
{
    namespace io
    {
        // Directory functions

        void create_work_directory(const char *folderpath);



        // Printing functions

        template <class T>
        void print_vector(std::vector<T> &vec);

        template <class T1, class T2>
        void print_vector_pairs(std::vector<std::pair<T1, T2>> &vec);

        template <class T>
        void print_matrix(std::vector<std::vector<T>> &mat);

        template <class T1, class T2>
        void print_matrix_pairs(std::vector<std::vector<std::pair<T1, T2>>> &mat);



        // Reading functions

        bool read_state(
            std::string file_path_name,
            std::vector<std::vector<double>> &X
        );

        std::list<std::pair<int, int>> read_total_task_time_pairs(
            std::string file_path_name
        );

        std::list<std::pair<std::pair<int, int>, double>> read_time_pairs_resultsfile(
            std::string file_path_name
        );

        bool get_total_ready_results(
            std::list<std::pair<std::pair<int, int>, double>>& total_ready_results,
            std::string output_folder_path
        );



        // Writing functions

        bool write_matrix_result(
            std::string output_file_path_name_matrix,
            std::vector<std::vector<double>>& mat
        );

        bool write_list_result(
            std::string output_file_path_name_list,
            std::pair<std::pair<int, int>, double>& result_pair
        );

        bool write_matching_result(
            std::string output_file_path_name_matching,
            std::vector<std::vector<double>>& Xa_vertex,
            std::vector<std::vector<double>>& Xb_vertex,
            std::vector<std::vector<double>>& cost_adjmatrix,
            double distance,
            std::vector<int>& q
        );
    }
};
