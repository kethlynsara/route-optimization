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
    double lat1 = to_radians(a.x), lon1 = to_radians(a.y);
    double lat2 = to_radians(b.x), lon2 = to_radians(a.y);
    double dlat = lat2 - lat1, dlon = lon2 - lon1;

    double a_hav = pow(sin(dlat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(dlon / 2), 2);
    double c = 2 * atan2(sqrt(a_hav), sqrt(1 - a_hav));

    return EARTH_RADIUS * c; // Distância em km
}

// Função para ler os dados da instância
void read_data(vector<Node> &nodes, int &dimension, string &edge_weight_type, string &output_file)
{
    cout << "Digite o nome do arquivo de saída: ";
    cin >> output_file;

    cout << "Forneça os dados da instância...\n";
    string line;
    while (getline(cin, line))
    {
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        if (line.rfind("DIMENSION", 0) == 0)
        {
            dimension = stoi(line.substr(line.find(":") + 1));
            cout << "Dimensão da instância: " << dimension << " pontos\n";
        }
        else if (line.rfind("EDGE_WEIGHT_TYPE", 0) == 0)
        {
            edge_weight_type = line.substr(line.find(":") + 1);
            edge_weight_type.erase(0, edge_weight_type.find_first_not_of(" \t"));
            cout << "Tipo de distância: " << edge_weight_type << "\n";
        }
        else if (line == "NODE_COORD_SECTION")
        {
            cout << "Lendo coordenadas dos nós...\n";
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

// Função para criar matriz de adjacência
vector<vector<double>> compute_adjacency_matrix(const vector<Node> &nodes, const string &edge_weight_type)
{
    cout << "Calculando matriz de adjacência...\n";
    int n = nodes.size();
    vector<vector<double>> adjacency_matrix(n, vector<double>(n, 0));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            adjacency_matrix[i][j] = (edge_weight_type == "EUC_2D") ? euclidean_distance(nodes[i], nodes[j]) : haversine_distance(nodes[i], nodes[j]);
        }
    }

    cout << "Matriz de adjacência calculada com sucesso!\n";
    return adjacency_matrix;
}

// Algoritmo (Nearest Neighbor)
vector<int> nearest_neighbor(const vector<vector<double>> &adj_matrix)
{
    cout << "Executando o algoritmo (Nearest Neighbor)...\n";
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

    cout << "Algoritmo guloso concluído!\n";
    return tour;
}

// Função para calcular o custo do percurso
double compute_tour_cost(const vector<int> &tour, const vector<vector<double>> &adj_matrix)
{
    double cost = 0;
    int n = tour.size();
    for (int i = 0; i < n; ++i)
    {
        cost += adj_matrix[tour[i]][tour[(i + 1) % n]];
    }
    return cost;
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
                continue;

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

int main()
{
    vector<Node> nodes;
    int dimension;
    string edge_weight_type, output_file;

    read_data(nodes, dimension, edge_weight_type, output_file);
    vector<vector<double>> adj_matrix = compute_adjacency_matrix(nodes, edge_weight_type);

    auto start = high_resolution_clock::now();
    vector<int> tour = nearest_neighbor(adj_matrix);
    double initial_cost = compute_tour_cost(tour, adj_matrix);

    cout << "Iniciando busca local (2-opt)...\n";
    while (two_opt(tour, adj_matrix));
    cout << "Busca local concluída!\n";

    double final_cost = compute_tour_cost(tour, adj_matrix);
    auto stop = high_resolution_clock::now();

    double percentage_improvement = 100.0 * (initial_cost - final_cost) / initial_cost;
    auto duration = duration_cast<milliseconds>(stop - start).count();

    // **Saída formatada corretamente no arquivo**
    cout << "Gravando solução no arquivo " << output_file << "...\n";
    ofstream output(output_file);
    for (size_t i = 0; i < tour.size(); i++)
    {
        output << "v_" << nodes[tour[i]].id;
        if (i < tour.size() - 1)
            output << " ";
    }
    output.close();
    cout << "Solução salva!\n";

    // **Resumo final**
    cout << "-------------------------------\n";
    cout << "Solução inicial: " << initial_cost << "\n";
    cout << "Solução final: " << final_cost << "\n";
    cout << "Melhoria: " << percentage_improvement << "%\n";
    cout << "Tempo total: " << duration << " ms\n";
    cout << "-------------------------------\n";

    return 0;
}
