#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>

using namespace std;
using namespace std::chrono;

const double EARTH_RADIUS = 6371.0;

struct Node
{
    int id;
    double x, y;
};

double euclidean_distance(const Node &a, const Node &b)
{
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

double to_radians(double degrees)
{
    return degrees * M_PI / 180.0;
}

double haversine_distance(const Node &a, const Node &b)
{
    double lat1 = to_radians(a.x), lon1 = to_radians(a.y);
    double lat2 = to_radians(b.x), lon2 = to_radians(b.y);

    double dlat = lat2 - lat1, dlon = lon2 - lon1;
    double a_hav = pow(sin(dlat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(dlon / 2), 2);
    double c = 2 * atan2(sqrt(a_hav), sqrt(1 - a_hav));

    return EARTH_RADIUS * c;
}

void read_data(vector<Node> &nodes, int &dimension, string &edge_weight_type, string &output_file)
{
    cout << "Digite o nome do arquivo de saída: " << endl;
    cin >> output_file;
    string line;
    while (getline(cin, line))
    {
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        if (line.rfind("DIMENSION", 0) == 0)
        {
            dimension = stoi(line.substr(line.find(":") + 1));
        }
        else if (line.rfind("EDGE_WEIGHT_TYPE", 0) == 0)
        {
            edge_weight_type = line.substr(line.find(":") + 1);
            edge_weight_type.erase(0, edge_weight_type.find_first_not_of(" \t"));
        }
        else if (line == "NODE_COORD_SECTION")
        {
            for (int i = 0; i < dimension; ++i)
            {
                Node node;
                cin >> node.id >> node.x >> node.y;
                nodes.push_back(node);
            }
        }
        else if (line == "EOF")
        {
            break;
        }
    }
}

vector<vector<double>> compute_adjacency_matrix(const vector<Node> &nodes, const string &edge_weight_type)
{
    int n = nodes.size();
    vector<vector<double>> adj_matrix(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            adj_matrix[i][j] = (edge_weight_type == "EUC_2D") ? euclidean_distance(nodes[i], nodes[j]) : haversine_distance(nodes[i], nodes[j]);
        }
    }
    return adj_matrix;
}

vector<int> nearest_neighbor(const vector<vector<double>> &adj_matrix)
{
    int n = adj_matrix.size();
    vector<int> tour, visited(n, false);
    int current = 0;
    tour.push_back(current);
    visited[current] = true;
    for (int i = 1; i < n; ++i)
    {
        int next = -1;
        double min_dist = numeric_limits<double>::max();
        for (int j = 0; j < n; ++j)
        {
            if (!visited[j] && adj_matrix[current][j] < min_dist)
            {
                min_dist = adj_matrix[current][j];
                next = j;
            }
        }
        tour.push_back(next);
        visited[next] = true;
        current = next;
    }
    return tour;
}

double max_edge_in_tour(const vector<int> &tour, const vector<vector<double>> &adj_matrix)
{
    double max_edge = 0;
    int n = tour.size();
    for (int i = 0; i < n; ++i)
    {
        max_edge = max(max_edge, adj_matrix[tour[i]][tour[(i + 1) % n]]);
    }
    return max_edge;
}

bool optimize_tour(vector<int> &tour, const vector<vector<double>> &adj_matrix)
{
    int n = tour.size();
    bool improved = false;
    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = i + 2; j < n; ++j)
        {
            double d1 = adj_matrix[tour[i]][tour[i + 1]];
            double d2 = adj_matrix[tour[j]][tour[(j + 1) % n]];
            double new_d1 = adj_matrix[tour[i]][tour[j]];
            double new_d2 = adj_matrix[tour[i + 1]][tour[(j + 1) % n]];
            if (max(new_d1, new_d2) < max(d1, d2))
            {
                reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                improved = true;
            }
        }
    }
    return improved;
}

void solve_instance(vector<Node> &nodes, const string &edge_weight_type, const string &output_file)
{
    vector<vector<double>> adj_matrix = compute_adjacency_matrix(nodes, edge_weight_type);
    auto start = high_resolution_clock::now();
    vector<int> tour = nearest_neighbor(adj_matrix);
    double initial_max_edge = max_edge_in_tour(tour, adj_matrix);
    while (optimize_tour(tour, adj_matrix))
        ;
    double final_max_edge = max_edge_in_tour(tour, adj_matrix);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start).count();
    ofstream output(output_file);
    for (size_t i = 0; i < tour.size(); i++)
    {
        output << "v_" << nodes[tour[i]].id << (i < tour.size() - 1 ? " " : "");
    }
    output.close();
    cout << "-------------------------------\n";
    cout << "Solução inicial: " << initial_max_edge << "\n";
    cout << "Solução final: " << final_max_edge << "\n";
    cout << "Desvio percentual: " << (100 * (initial_max_edge - final_max_edge) / initial_max_edge) << "%\n";
    cout << "Tempo total: " << duration << " ms\n";
    cout << "-------------------------------\n";
}

int main()
{
    string output_file;
    vector<Node> nodes;
    int dimension;
    string edge_weight_type;
    read_data(nodes, dimension, edge_weight_type, output_file);
    solve_instance(nodes, edge_weight_type, output_file);
    return 0;
}
