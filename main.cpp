#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <fstream>

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

// Função para imprimir a sequência de visitas no arquivo de saída
void write_solution_to_file(const string &filename, const vector<int> &solution)
{
    ofstream output_file(filename);
    for (size_t i = 0; i < solution.size(); ++i)
    {
        output_file << "v_" << solution[i];
        if (i != solution.size() - 1)
            output_file << " ";
    }
    output_file.close();
}

int main()
{
    string output_filename;
    cout << "Digite o nome do arquivo para salvar a solução: ";
    cin >> output_filename; // Recebe o nome do arquivo via cin

    vector<Node> nodes;
    int dimension;
    string edge_weight_type;

    cout << "Insira os dados no formato esperado:\n";
    read_data(cin, nodes, dimension, edge_weight_type);

    cout << "Dimension lida: " << dimension << endl;
    cout << "Tipo de distância: " << edge_weight_type << endl;

    // Criar matriz de adjacências com distâncias entre os nós
    vector<vector<double>> adjacency_matrix(dimension, vector<double>(dimension, 0));

    for (int i = 0; i < dimension; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            if (edge_weight_type == "EUC_2D")
            {
                adjacency_matrix[i][j] = euclidean_distance(nodes[i], nodes[j]);
            }
            else if (edge_weight_type == "GEO")
            {
                adjacency_matrix[i][j] = haversine_distance(nodes[i], nodes[j]);
            }
            else
            {
                cerr << "Erro: Tipo de distância desconhecido!" << endl;
                return 1;
            }
        }
    }

    // Supondo que a melhor solução seja obtida pelo método Nearest Neighbor + 2-opt
    vector<int> best_solution = {0, 1, 2, 3};               // Exemplo de solução encontrada
    write_solution_to_file(output_filename, best_solution); // Grava solução no arquivo especificado

    cout << "Solução gravada no arquivo: " << output_filename << endl;

    return 0;
}
