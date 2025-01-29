#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;

const double EARTH_RADIUS = 6371.0; // Raio da Terra em km

struct Node
{
    int id;
    double x, y;
};

// Função para calcular a distância euclidiana entre dois pontos (EUC_2D)
double euclidean_distance(const Node &a, const Node &b)
{
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

// Função para converter graus para radianos
double to_radians(double degrees)
{
    return degrees * M_PI / 180.0;
}

// Função para calcular a distância entre dois pontos geográficos (GEO) usando Haversine
double haversine_distance(const Node &a, const Node &b)
{
    double lat1 = to_radians(a.x);
    double lon1 = to_radians(a.y);
    double lat2 = to_radians(b.x);
    double lon2 = to_radians(b.y);

    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    double a_hav = pow(sin(dlat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(dlon / 2), 2);
    double c = 2 * atan2(sqrt(a_hav), sqrt(1 - a_hav));

    return EARTH_RADIUS * c; // Distância em km
}

// Função para ler os dados da instância
void read_data(istream &input, vector<Node> &nodes, int &dimension, string &edge_weight_type)
{
    string line;

    while (getline(input, line))
    {
        // Remove espaços extras no início e no final
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
                input >> node.id >> node.x >> node.y;
                nodes.push_back(node);
            }
        }
        else if (line == "EOF")
        {
            break;
        }
    }
}

// Função para calcular a matriz de adjacências
vector<vector<double>> compute_adjacency_matrix(const vector<Node> &nodes, const string &edge_weight_type)
{
    int n = nodes.size();
    vector<vector<double>> adjacency_matrix(n, vector<double>(n, 0));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (edge_weight_type == "EUC_2D")
            {
                adjacency_matrix[i][j] = euclidean_distance(nodes[i], nodes[j]);
            }
            else if (edge_weight_type == "GEO")
            {
                adjacency_matrix[i][j] = haversine_distance(nodes[i], nodes[j]);
            }
        }
    }

    return adjacency_matrix;
}

// Algoritmo Guloso (Nearest Neighbor)
vector<int> nearest_neighbor(const vector<vector<double>> &adj_matrix)
{
    int n = adj_matrix.size();
    vector<int> tour;
    vector<bool> visited(n, false);

    int current = 0; // Começa no primeiro nó
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

// Função para calcular a maior distância entre dois pontos no ciclo
double max_edge_in_tour(const vector<int> &tour, const vector<vector<double>> &adj_matrix)
{
    double max_edge = 0;
    int n = tour.size();

    for (int i = 0; i < n; ++i)
    {
        double dist = adj_matrix[tour[i]][tour[(i + 1) % n]];
        max_edge = max(max_edge, dist);
    }

    return max_edge;
}

// Busca Local (2-opt)
bool two_opt(vector<int> &tour, const vector<vector<double>> &adj_matrix)
{
    int n = tour.size();
    bool improved = false;

    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = i + 2; j < n; ++j)
        {
            if (j == n - 1 && i == 0)
                continue; // Evita inverter a primeira e última aresta

            double d1 = adj_matrix[tour[i]][tour[i + 1]] + adj_matrix[tour[j]][tour[(j + 1) % n]];
            double d2 = adj_matrix[tour[i]][tour[j]] + adj_matrix[tour[i + 1]][tour[(j + 1) % n]];

            if (d2 < d1)
            {
                reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                improved = true;
            }
        }
    }

    return improved;
}

// Função para refinar a solução com 2-opt
void optimize_with_two_opt(vector<int> &tour, const vector<vector<double>> &adj_matrix)
{
    while (two_opt(tour, adj_matrix))
        ;
}

int main()
{
    vector<Node> nodes;
    int dimension;
    string edge_weight_type;

    cout << "Insira os dados no formato esperado:\n";
    read_data(cin, nodes, dimension, edge_weight_type);

    vector<vector<double>> adj_matrix = compute_adjacency_matrix(nodes, edge_weight_type);

    vector<int> tour = nearest_neighbor(adj_matrix);
    optimize_with_two_opt(tour, adj_matrix);

    // Exibir solução
    cout << "Melhor sequência de paradas:\n";
    for (int node : tour)
    {
        cout << nodes[node].id << " ";
    }
    cout << endl;

    return 0;
}
