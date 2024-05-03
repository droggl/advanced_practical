#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <random>
#include <queue>
#include <chrono>
#include <stack>
#include <bits/stdc++.h>

using namespace std;

typedef unsigned int NodeID;
typedef unsigned int EdgeID;
typedef unsigned int PartitionID;

int LOG_STEPS = 10;

bool randomBoolean() {
  static default_random_engine generator(std::random_device{}());
  // With p = 0.5 you get equal probability for true and false
  static bernoulli_distribution distribution(0.5);
  return distribution(generator);
}

template <typename T>
T randomFrom(const T min, const T max)
{
    static std::random_device rdev;
    static std::default_random_engine re(rdev());
    typedef typename std::conditional<
        std::is_floating_point<T>::value,
        std::uniform_real_distribution<T>,
        std::uniform_int_distribution<T>>::type dist_type;
    dist_type uni(min, max);
    return static_cast<T>(uni(re));
}


struct Node {
    vector<NodeID> adjacents;
};

class Graph {

    public:
        int n, m;
        vector<Node> nodes;
        double gamma, alpha, v;

        Graph(){}

        int read_edge_list_from_file(string filename) {
            ifstream file(filename);

            if (!file.is_open()) {
                cerr << "Error: Could not open file " << filename << endl;
                return 1;
            }

            string line;

            // Read the header to get the number of nodes and edges
            while (getline(file, line)) {
                if (line.empty() || line[0] == '%') continue; // Ignore comments and empty lines
                stringstream ss(line);
                ss >> n >> n >> m;

                break;
            }
            nodes.resize(n);
            // Read the edges
            int new_m = 0;
            while(getline(file, line)) {
                int from_node_id, to_node_id;

                stringstream ss(line);
                ss >> from_node_id >> to_node_id;
                // Decrement IDs by 1
                --from_node_id;
                --to_node_id;

                if (find(nodes[from_node_id].adjacents.begin(), nodes[from_node_id].adjacents.end(), to_node_id) == nodes[from_node_id].adjacents.end()) {
                    nodes[from_node_id].adjacents.push_back(to_node_id);
                    nodes[to_node_id].adjacents.push_back(from_node_id);
                    new_m ++;
                }

            }
            m = new_m;

            return 0;
        }

        int read_graph_from_metis_file(string filename) {
            ifstream file(filename);
            if (!file.is_open()) {
                cerr << "Error: Could not open file " << filename << endl;
                return 1;
            }

            string line;
            // Skip comment lines
            while (getline(file, line) && line[0] == '%');

            // Read number of vertices and edges
            stringstream ss(line);
            ss >> n >> m;

            nodes.resize(n);

            int node_id = 0;
            while (getline(file, line)) {
                if (line[0] == '%')
                    continue;
                stringstream ss(line);
                int neighbor;
                while (ss >> neighbor) {
                    nodes[node_id].adjacents.push_back(neighbor-1);
                }
                node_id++;
            }

            file.close();


            return 0;
        }

        static Graph generate_hp_random_graph(int num_nodes, int num_partitions, double intra_prob, double inter_prob) {
            Graph graph;
            graph.n = num_nodes;
            graph.nodes.resize(graph.n);
            vector<int> partitions(graph.n);

            int partition = 0;
            for (int i=0; i<num_nodes; i++) {
                partitions[i] = partition;
                partition = (partition + 1) % num_partitions;
            }

            random_device rd;
            mt19937 gen(rd());
            uniform_real_distribution<> dis(0.0, 1.0);

            int m_new = 0;

            // Generation of edges according to probabilities
            for (int i = 0; i < num_nodes; ++i) {
                for (int j = i + 1; j < num_nodes; ++j) {
                    if (partitions[i] == partitions[j]) {
                        if (dis(gen) < intra_prob) {
                            graph.nodes[i].adjacents.push_back(j);
                            graph.nodes[j].adjacents.push_back(i);
                            m_new++;
                        }
                    } else {
                        if (dis(gen) < inter_prob) {
                            graph.nodes[i].adjacents.push_back(j);
                            graph.nodes[j].adjacents.push_back(i);
                            m_new++;
                        }
                    }
                }
            }
            graph.m = m_new;

            return graph;
        }
};



template<typename T> void printElement(T t, const int& width){
    const char separator = ' ';
    cout << left << setw(width) << setfill(separator) << setprecision(3) << t << "|";
}

const int n_width = 7;
const int m_width = 10;
const int k_width = 4;
const int lambda_width = 8;
const int rho_width = 7;

class Partition {
    public:
        Graph graph;
        vector<PartitionID> partition;
        vector<int> nodes_in_partition;
        int k;

        double v = 1.1;
        double gamma = 1.5;
        double alpha, lambda, rho;
        double cutsize;

        Partition(Graph g, int num_partitions) : graph(g), k(num_partitions), partition(graph.n, -1), nodes_in_partition(num_partitions, 0) {
            alpha = (double) graph.m * pow(k, 0.5) / pow(graph.n, gamma);
        };

        void stream_partition(string permutation, string status="progress", bool enforce_balance_constraint=true) {
            int progress = 0;
            if (status=="progress"){
                cout << "|" << flush;
            }

            int capacity = v * (graph.n/ ((float) k));
            vector<NodeID> node_ordering;
            get_node_ordering(node_ordering, permutation);


            for (int i = 0; i<graph.n; i++) {
                double max_gain = std::numeric_limits<float>::lowest();
                double gain;
                PartitionID best_partition = randomFrom(0, k-1);
                NodeID node_id = node_ordering[i];
                for(PartitionID p_id=0; p_id<k; p_id++) {
                    // Check if new node in partition would exceed threshold v
                    if (enforce_balance_constraint && nodes_in_partition[p_id] >= capacity) {
                        continue;
                    }
                    gain = fennel_gain(node_id, p_id, nodes_in_partition[p_id]);
                    if (gain > max_gain || (randomBoolean() && gain == max_gain)) {
                        max_gain = gain;
                        best_partition = p_id;
                    }
                }
                partition[node_id] = best_partition;
                nodes_in_partition[best_partition]++;

                if (status == "progress") {
                    log_progress(i, progress);
                }
            }

            if (status == "progress") {
                cout << "|\n" << endl;
            }

            eval_partition();
        }

        void print_partition_evaluation() {
            Partition::print_partition_evaluation_title();
            Partition::print_partition_evaluation(graph.n, graph.m, k, cutsize, 100*lambda, rho);
        }

        static void print_partition_evaluation_title() {
            printElement("n", n_width);
            printElement("m", m_width);
            printElement("k", k_width);
            printElement("cutsize", n_width);
            printElement("lambda", lambda_width);
            printElement("rho", rho_width);
            cout << endl << "---------------------------------------------------" << endl;
        }

        static void print_partition_evaluation(int n, int m, int num_partition, int cutsize, double lambda, double rho) {
            printElement(n, n_width);
            printElement(m, m_width);
            printElement(num_partition, k_width);
            printElement(cutsize, n_width);
            printElement(lambda, lambda_width);
            printElement(rho, rho_width);
            cout << endl;
        }

    private:

        void get_node_ordering(vector<NodeID>& node_ordering, string permutation) {
            if (permutation == "standard" || permutation == "random") {
                node_ordering.resize(graph.n);
                for (NodeID i = 0; i < graph.n; ++i) {
                    node_ordering[i] = i;
                }
                if (permutation == "random") {
                    random_shuffle(node_ordering.begin(), node_ordering.end());
                }
            } else if (permutation == "dfs" || permutation == "bfs") {
                random_device rd;
                mt19937 gen(rd());
                uniform_int_distribution<> dis(0, graph.n-1);
                int start_node = dis(gen);
                if (permutation == "dfs") {
                    dfs(node_ordering, start_node);
                } else {
                    bfs(node_ordering, start_node);
                }
            }
        }

        void dfs(vector<NodeID>& ordering, NodeID start_node_id) {
            vector<bool> already_visited(graph.n, false);
            stack<NodeID> stack;
            stack.push(start_node_id);

            while (!stack.empty()) {
                NodeID current_node_id = stack.top();
                stack.pop();

                if (!already_visited[current_node_id]) {
                    already_visited[current_node_id] = true;
                    ordering.push_back(current_node_id);

                    for (NodeID adjacent_id : graph.nodes[current_node_id].adjacents) {
                        if (!already_visited[adjacent_id]) {
                            stack.push(adjacent_id);
                        }
                    }
                }
            }

            NodeID node_id=0;
            for (bool visited_node : already_visited) {
                if (!visited_node) {
                    ordering.push_back(node_id);
                }
                node_id++;
            }

        }

        void bfs(vector<NodeID>& ordering, NodeID start_node_id) {
            queue<NodeID> queue;
            vector<bool> visited(graph.nodes.size(), false);

            queue.push(start_node_id);
            visited[start_node_id] = true;
            while (!queue.empty()) {
                NodeID current_node_id = queue.front();
                queue.pop();

                ordering.push_back(current_node_id);
                for (NodeID adjacent_id : graph.nodes[current_node_id].adjacents) {
                    visited[adjacent_id];
                    if (!visited[adjacent_id]) {
                        queue.push(adjacent_id);
                        visited[adjacent_id] = true;
                    }
                }
            }

            NodeID node_id=0;
            for (bool visited_node : visited) {
                if (!visited_node) {
                    ordering.push_back(node_id);
                }
                node_id++;
            }
        }

        double fennel_gain(NodeID node_id, PartitionID p_id, int nodes_in_partition) {
            int neighbours_in_partition = 0;
            for (NodeID adj: graph.nodes[node_id].adjacents) {
                if (partition[adj] == p_id) {
                    neighbours_in_partition++;
                }
            }

            double intra_partition_cost = alpha * gamma * sqrt(nodes_in_partition);
            return neighbours_in_partition - intra_partition_cost;
        }

        void compute_cut_size() {
            cutsize = 0;
            for (NodeID node_id = 0; node_id < graph.n; node_id++) {
                for (NodeID neighbor_id : graph.nodes[node_id].adjacents) {
                    if (partition[neighbor_id] != partition[node_id]) {
                        cutsize++;
                    }
                }
            }
            cutsize /= 2; // Divide by 2 to account for doubled edges
        }

        void eval_partition() {

            compute_cut_size();
            lambda = cutsize / graph.m;
            int max_load = 0;
            for (PartitionID i=0; i<k; i++) {
                int load = count(partition.begin(), partition.end(), i);
                if (load > max_load) {
                    max_load = load;
                }
            }
            rho = max_load / (partition.size() / (double) k);
        }

        void log_progress(int i, int& progress) {
            if ((i % (int)(1.0/LOG_STEPS * graph.n)) == 0) {
                cout << "o" << flush;
                progress++;
            }
        }
};


void evaluate_hp_random_graphs(string ordering) {
    int num_nodes = 5000;
    double intra_prob = 0.8;
    double inter_prob = 0.5;

    int rounds = 5;
    vector<int> v_k = {4, 8, 16, 32, 64, 128};
    double avg_lambda, avg_rho, avg_m, avg_cutsize;
    cout << "Ordering: " << ordering << endl;
    Partition::print_partition_evaluation_title();

    for (int num_partitions : v_k) {
        avg_lambda = avg_rho = avg_m = avg_cutsize = 0;

        for (int i=0; i<rounds; i++) {
            Graph random_graph = Graph::generate_hp_random_graph(num_nodes, num_partitions, intra_prob, inter_prob);
            Partition p = Partition(random_graph, num_partitions);
            p.stream_partition(ordering, "none", false);

            avg_lambda += p.lambda;
            avg_rho += p.rho;
            avg_m += random_graph.m;
            avg_cutsize += p.cutsize;
        }

        avg_lambda =  100 * avg_lambda / rounds;
        avg_rho /= rounds;
        avg_m /= rounds;
        avg_cutsize /=rounds;

        Partition::print_partition_evaluation(num_nodes, (int) avg_m, num_partitions, avg_cutsize, avg_lambda, avg_rho);
    }
}

bool ORDERING_STANDARD = true;
bool ORDERING_RANDOM = true;
bool ORDERING_DFS = true;
bool ORDERING_BFS = true;

int main(int argc, char* argv[]) {


    if (argc == 1) {
        evaluate_hp_random_graphs("random");
    } else if (argc != 3) {
            cerr << "Usage: " << argv[0] << " <filename>" << "(optional) <k>" << endl;
            return 1;
    } else if(argc == 3) {
        int k =  atoi(argv[2]);

        string filename = argv[1];


        Graph graph = Graph();
        int ret_val = graph.read_graph_from_metis_file(filename);
        if (ret_val != 0) {
            return ret_val;
        }

        string ordering;

        if (ORDERING_STANDARD){
            ordering = "standard";
            cout << endl << "Ordering: " << ordering << endl;
            Partition p_std = Partition(graph, k);
            p_std.stream_partition(ordering, "none");
            p_std.print_partition_evaluation();
        } if (ORDERING_RANDOM) {
            ordering = "random";
            cout << endl << "Ordering: " << ordering << endl;
            Partition p_random = Partition(graph, k);
            p_random.stream_partition(ordering, "none");
            p_random.print_partition_evaluation();
        } if (ORDERING_DFS) {
            ordering = "dfs";
            cout << endl << "Ordering: " << ordering << endl;
            Partition p_dfs = Partition(graph, k);
            p_dfs.stream_partition(ordering, "none");
            p_dfs.print_partition_evaluation();
        } if (ORDERING_BFS) {
            ordering = "bfs";
            cout << endl << "Ordering: " << ordering << endl;
            Partition p_bfs = Partition(graph, k);
            p_bfs.stream_partition(ordering, "none");
            p_bfs.print_partition_evaluation();
        }

        cout << endl;
    }

    return 0;
}
