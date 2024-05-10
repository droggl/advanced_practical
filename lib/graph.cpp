#include "graph.h"

int Graph::read_edge_list_from_file(string filename) {
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
    while (getline(file, line)) {
        int from_node_id, to_node_id;

        stringstream ss(line);
        ss >> from_node_id >> to_node_id;
        // Decrement IDs by 1
        --from_node_id;
        --to_node_id;

        if (find(nodes[from_node_id].adjacents.begin(), nodes[from_node_id].adjacents.end(), to_node_id) == nodes[from_node_id].adjacents.end()) {
            nodes[from_node_id].adjacents.push_back(to_node_id);
            nodes[to_node_id].adjacents.push_back(from_node_id);
            new_m++;
        }
    }
    m = new_m;

    return 0;
}

int Graph::read_graph_from_metis_file(string filename) {
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
            nodes[node_id].adjacents.push_back(neighbor - 1);
        }
        node_id++;
    }

    file.close();

    return 0;
}

Graph Graph::extract_subgraph(vector<NodeID>& ordering, NodeID start_node, int buffer_size, vector<NodeID>& local_to_global) const { //,unordered_map<NodeID, NodeID>& local_to_global) const {
    Graph subgraph;
    // Set subgraph.n to buffer_size if does not extend n
    subgraph.n = start_node + buffer_size > n ? n - start_node : buffer_size;
    subgraph.nodes.resize(subgraph.n);
    local_to_global.resize(subgraph.n);
    subgraph.m = 0;

    if (buffer_size == 1) {
        local_to_global[0] = start_node;
        return subgraph;
    }

    unordered_map<NodeID, NodeID> global_to_local;

    // Create global to local (and local to global) assignment
    for (NodeID i = 0; i < subgraph.n; ++i) {
        NodeID global_id = ordering[start_node+i];
        global_to_local[global_id] = i;
        local_to_global[i] = global_id;
    }

    // Extract adjacency list of the node in the subgraph
    for (NodeID i = 0; i < subgraph.n; ++i) {
        NodeID global_node_id = ordering[start_node+i];
        const Node& global_node = nodes[global_node_id];
        Node& local_node = subgraph.nodes[i];

        for (NodeID global_adj : global_node.adjacents) {
            // Check if neighbour is part of subgraph
            if (global_to_local.find(global_adj) != global_to_local.end()) {
                local_node.adjacents.push_back(global_to_local[global_adj]);
                subgraph.m++;
            }
        }
    }

    subgraph.m /= 2;
    return subgraph;
}


Graph Graph::generate_hp_random_graph(int num_nodes, int num_partitions, double intra_prob, double inter_prob) {
    Graph graph;
    graph.n = num_nodes;
    graph.nodes.resize(graph.n);
    vector<int> partitions(graph.n);

    int partition = 0;
    for (int i = 0; i < num_nodes; i++) {
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