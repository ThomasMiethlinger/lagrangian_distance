// Author:     Thomas Miethlinger B.Sc.
// Date:       01.2020
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include "../include/io.hpp"

#include <iostream> // std::cout, std::cerr
#include <fstream> // std::ifstream, std::ofstream
#include <list> // std::list
#include <string> // std::string
#include <vector> // std::vector
#include <utility> // std::pair

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h> //dirent, DIR

// Boost
#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::list;
using std::pair;
using std::string;
using std::vector;

typedef std::vector<int> VI;
typedef std::vector<std::vector<double>> VVD;

// Directory functions

// Create directory or write error unless the directory already exists
void util::io::create_work_directory(const char *folderpath)
{
    if(mkdir(folderpath, 0777) == -1 && errno != 17) // errno == 17 -> directory already exists
        std::cerr << "Error: " << strerror(errno) << ". " << folderpath << std::endl;
}



// Printing functions

template <class T>
void util::io::print_vector(vector<T> &vec)
{
    for(int i = 0; i < vec.size(); i++)
    {
        cout << vec[i] << endl;
    }
    cout << endl;
}

template <class T1, class T2>
void util::io::print_vector_pairs(vector<pair<T1, T2>> &vec)
{    
    for(int i = 0; i < vec.size(); i++)
    {
        cout << vec[i].first << " " << vec[i].second << endl;
    }
    cout << endl;
}

template <class T>
void util::io::print_matrix(vector<vector<T>> &mat)
{
    for(int i = 0; i < mat.size(); i++)
    {
        for(int j = 0; j < mat[i].size(); j++)
        {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

template <class T1, class T2>
void util::io::print_matrix_pairs(vector<vector<pair<T1, T2>>> &mat)
{    
    for(int i = 0; i < mat.size(); i++)
    {
        for(int j = 0; j < mat[i].size(); j++)
        {
            cout << mat[i][j].first << " " << mat[i][j].second << "   ";
        }
        cout << endl;
    }
    cout << endl;
}



// Reading functions

// Designed for a Lammps/Liggghts dump file of the following form:
// id x y z vx vy vz
// X is a vector of size N, consisting of vectors of size 6. Thus, X is a N times 6 double matrix
bool util::io::read_state(string file_path_name, vector<vector<double>> &X)
{
    int N = X.size();

    bool result = true;

    std::ifstream filestream;
    filestream.open(file_path_name);
    if (filestream.is_open())
    {
        string line;
        vector<string> parts;
        for (int i = -9; std::getline(filestream, line) && i < N; i++)
        {
            if (i >= 0)
            {
                boost::split(parts, line, boost::is_any_of(" "));
                vector<double> x(6);
                for (int j = 0; j < 6; j++)
                {
                    x[j] = std::stod(parts[j + 1]);
                }
                X[i] = x;
            }
        }
        filestream.close();
    }
    else
    {
        cout << "Error: Could not read file " << file_path_name << "!" << endl;
        result = false;
    }

    return result;
}

// Read time-pairs defined by existing result files
list<pair<int, int>> util::io::read_total_task_time_pairs(string file_path_name)
{
    list<pair<int, int>> total_task_time_pairs;
    pair<int, int> time_pair;
    int t1, t2;
    std::ifstream filestream;
    filestream.open(file_path_name);
    if(filestream.is_open())
    {
        string line;
        vector<string> parts;
        for(int i = 0; std::getline(filestream, line); i++)
        {
            boost::split(parts, line, boost::is_any_of(" "));
            t1 = std::stoi(parts[0]);
            t2 = std::stoi(parts[1]);
            time_pair = std::make_pair(t1, t2);
            total_task_time_pairs.push_back(time_pair);
        }
        filestream.close();
    }
    else
    {
        cout << "Error: Could not read file " << file_path_name << "!" << endl;
    }

    return total_task_time_pairs;
}

// Read time-pairs defined by existing result files
list<pair<pair<int, int>, double>> util::io::read_time_pairs_resultsfile(string file_path_name)
{
    list<pair<pair<int, int>, double>> inputlist_timepairs;
    pair<pair<int, int>, double> row;
    pair<int, int> time_pair;
    int t1, t2;
    double d;
    std::ifstream filestream;
    filestream.open(file_path_name);
    if(filestream.is_open())
    {
        string line;
        vector<string> parts;
        for(int i = 0; std::getline(filestream, line); i++)
        {
            boost::split(parts, line, boost::is_any_of(" "));
            t1 = std::stoi(parts[0]);
            t2 = std::stoi(parts[1]);
            d  = std::stod(parts[2]);
            time_pair = std::make_pair(t1, t2);
            row = std::make_pair(time_pair, d);
            inputlist_timepairs.push_back(row);
        }
        filestream.close();
    }
    else
    {
        cout << "Error: Could not read file " << file_path_name << "!" << endl;
    }

    return inputlist_timepairs;
}

// Read all time-pairs defined by output folder path
bool util::io::get_total_ready_results(list<pair<pair<int, int>, double>>& total_ready_results, string output_folder_path)
{
    list<pair<pair<int, int>, double>> intermediate_result;

    struct dirent *de;
    DIR *dr = opendir(output_folder_path.c_str()); 

    if(dr == NULL)
    { 
        cout << "Could not open output directory " << output_folder_path << "!" << endl;
        return false;
    } 

    while((de = readdir(dr)) != NULL)
    {
        string file_path_name = output_folder_path + de->d_name;
        intermediate_result = read_time_pairs_resultsfile(file_path_name);
        total_ready_results.splice(total_ready_results.end(), intermediate_result);
    }
    closedir(dr);

    return true;
}



// Writing functions

bool util::io::write_matrix_result(std::string output_file_path_name_matrix, VVD& matrix)
{
    bool success = false;

    std::ofstream filestream;
    filestream.open(output_file_path_name_matrix);
    
    if(filestream.is_open())
    {
        for(int i = 0; i < matrix.size(); i++)
        {
            for(int j = 0; j < matrix[i].size(); j++)
            {
                filestream << matrix[i][j];
                if(j < matrix[i].size() - 1)
                    filestream << " ";
            }
            filestream << endl;
        }

        success = true;
        filestream.close();
    }
    else
    {
        success = false;
        cout << "Error: Could not write to file " << output_file_path_name_matrix << "!" << endl;
    }

    return success;
}

bool util::io::write_list_result(string output_file_path_name_list, pair<pair<int, int>, double>& result_pair)
{
    bool success = false;

    std::ofstream filestream;
    filestream.open(output_file_path_name_list, std::ios_base::app);
    
    if(filestream.is_open())
    {
        filestream << result_pair.first.first << " " << result_pair.first.second << " " << result_pair.second << endl;

        success = true;
        filestream.close();
    }
    else
    {
        success = false;
        cout << "Error: Could not write to file " << output_file_path_name_list << "!" << endl;
    }

    return success;
}

bool util::io::write_matching_result(string output_file_path_name_matching, VVD& Xa_vertex, VVD& Xb_vertex, VVD& cost_adjmatrix, double distance, VI& q)
{
    int N = Xa_vertex.size();
    int n = Xa_vertex[0].size();

    bool success = false;

    std::ofstream filestream;
    filestream.open(output_file_path_name_matching);
    
    if(filestream.is_open())
    {
        filestream << distance << std::endl;
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < n; j++)
            {
                filestream << Xa_vertex[i][j] << " ";
            }

            for(int j = 0; j < n; j++)
            {
                filestream << Xb_vertex[q[i]][j] << " ";
            }

            filestream << cost_adjmatrix[i][q[i]] << std::endl;
        }

        success = true;
        filestream.close();
    }
    else
    {
        success = false;
        cout << "Error: Could not write to file " << output_file_path_name_matching << "!" << endl;
    }

    return success;
}
