#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <queue>
#include <chrono>
#include <stack>
#include <bits/stdc++.h>

using namespace std;

typedef unsigned int NodeID;
typedef unsigned int EdgeID;
typedef unsigned int PartitionID;


vector<NodeID> random_permutation(int n) {
    vector<NodeID> permutation(n);
    for (int i = 0; i < n; ++i) {
        permutation[i] = i;
    }
    random_device rd;
    mt19937 g(rd());
    shuffle(permutation.begin(), permutation.end(), g);
    return permutation;
}



int LOG_STEPS = 10;

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
            while(getline(file, line)) {
                int from_node_id, to_node_id;

                stringstream ss(line);
                ss >> from_node_id >> to_node_id;
                // Decrement IDs by 1
                --from_node_id;
                --to_node_id;

                nodes[from_node_id].adjacents.push_back(to_node_id);
                nodes[to_node_id].adjacents.push_back(from_node_id);
            }

            return 0;
        }

    void traverse_iterative_dfs(vector<NodeID>& dfs_ordering, vector<bool>& already_visited, NodeID node_id) {
    stack<NodeID> stack;
    stack.push(node_id);

    while (!stack.empty()) {
        NodeID current_node_id = stack.top();
        stack.pop();

        if (!already_visited[current_node_id]) {
            already_visited[current_node_id] = true;
            dfs_ordering.push_back(current_node_id);

            for (NodeID adjacent_id : nodes[current_node_id].adjacents) {
                if (!already_visited[adjacent_id]) {
                    stack.push(adjacent_id);
                }
            }
        }
    }
}

    void log_progress(int i, int& progress) {
        // cout << "|";
        if ((i % (int)(1.0/LOG_STEPS * n)) == 0) {
            cout << "o" << flush;
            // for (int j=0; j<LOG_STEPS; j++) {
            //     string output = (j<progress) ? "o" : " ";
            //     cout << output;
            // }
            progress++;
        }
        // cout << "|" << endl;
    }

    void stream_partition(int k, string permutation, string heuristic="Fennel") {

        vector<PartitionID> partitions_fennel(n, -1);
        vector<int> nodes_in_partition(k, 0);

        v = 1.1;
        gamma = 1.5;
        alpha = sqrt(k) * m / pow(n, gamma);


        cout << "|" << flush;
        int progress = 0;

        vector<NodeID> node_order;
        if (permutation == "standard") {
            node_order.resize(n);
            for (NodeID i = 0; i < n; ++i) {
                node_order[i] = i;
            }
        } else if (permutation == "random"){
            node_order = random_permutation(n);
        } else if (permutation == "dfs") {
            vector<bool> already_visited(n, false);
            vector<NodeID> new_node_id(n);
            traverse_iterative_dfs(node_order, already_visited, 0);
        }

        int i = 0;
        for (NodeID node_id : node_order) {
            greedy_vertex_assignment(node_id, k, partitions_fennel, nodes_in_partition, heuristic);
            log_progress(i, progress);
            ++i;
        }

        cout << "|\n" << endl;

        print_partition_evaluation(partitions_fennel, k, heuristic, permutation);
    }

private:

    void greedy_vertex_assignment(NodeID node_id, int k, vector<PartitionID>& partitions, vector<int>& nodes_in_partition, string heuristic="Fennel") {
        int max_gain = -0xFFFFFFF;
        int neighbours_in_partition, gain, delta_c;
        PartitionID best_partition;
        bool partition_allowed;
        for(PartitionID p_id=0; p_id<k; p_id++) {
            // Check if new node in partition would exceed threshold v
            partition_allowed = nodes_in_partition[p_id] + 1 <= v * (n/k);
            if (partition_allowed) {
                 if (heuristic == "Fennel") {
                    gain = fennel_gain(partitions, node_id, p_id, nodes_in_partition[p_id]);
                } else if (heuristic == "LDG") {
                    gain = ldg_gain(partitions, node_id, p_id, k, nodes_in_partition[p_id]);
                }
                if (gain > max_gain) {
                    max_gain = gain;
                    best_partition = p_id;
                }
            }
        }

        partitions[node_id] = best_partition;
        nodes_in_partition[best_partition]++;
    }
    double fennel_gain(vector<PartitionID>& partitions, NodeID node_id, PartitionID p_id, int nodes_in_partition) {
        int neighbours_in_partition = 0;
        for (NodeID adj: nodes[node_id].adjacents) {
            if (partitions[adj] == p_id) {
                neighbours_in_partition++;
            }
        }

        double cost = alpha * gamma * pow(nodes_in_partition, gamma-1);
        return neighbours_in_partition - cost;
    }

    double ldg_gain(vector<PartitionID>& partitions, NodeID node_id, PartitionID p_id, int k, int nodes_in_partition) {
        int neighbours_in_partition = 0;
        for (NodeID adj: nodes[node_id].adjacents) {
            if (partitions[adj] == p_id) {
                neighbours_in_partition++;
            }
        }
        return (double) neighbours_in_partition * (1 - (nodes_in_partition / ((double) n/ (double) k)));
    }

    int compute_cut_size(const vector<PartitionID>& partitions) {
        int cut_size = 0;
        for (NodeID i = 0; i < n; ++i) {
            int cur_partition = partitions[i];
            for (NodeID neighbor : nodes[i].adjacents) {
                if (partitions[neighbor] != cur_partition) {
                    cut_size++;
                }
            }
        }
        return cut_size / 2; // Divide by 2 to account for doubled edges
    }

    void print_partition_evaluation(vector<PartitionID>& partitions, int k, string title, string permutation) {
        // Compute cut size
        int cut_size = compute_cut_size(partitions);

        // Compute Lambda
        double lambda = cut_size / (double) m;

        // Compute Rho
        int max_load = 0;
        for (PartitionID i=0; i<k; i++) {
            int load = count(partitions.begin(), partitions.end(), i);
            if (load > max_load) {
                max_load = load;
            }
        }
        double rho = max_load / (partitions.size() / (double) k);

        cout << "#-- " << title << ", ordering: " << permutation << " --#" << endl;
        cout << "N: " << n << ", M: " << m << ", k: " << k << endl;
        cout << "Cut Size: " << cut_size << endl;
        cout << "Lamda: " << lambda << endl;
        cout << "rho: "<< rho << endl << endl;
    }
};





int main(int argc, char* argv[]) {
    int k = 2;
    if (argc != 2 && argc != 3) {
        cerr << "Usage: " << argv[0] << " <filename>" << "(optional) <k>" << endl;
        return 1;
    } else if(argc == 3) {
        k =  atoi(argv[2]);
    }

    string filename = argv[1];



    Graph graph = Graph();
    int ret_val = graph.read_edge_list_from_file(filename);
    if (ret_val != 0) {
        return ret_val;
    }

    graph.stream_partition(k, "standard", "Fennel");
    graph.stream_partition(k, "random", "Fennel");
    graph.stream_partition(k, "dfs", "Fennel");

    return 0;
}
