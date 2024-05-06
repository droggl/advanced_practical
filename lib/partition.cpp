#include "partition.h"
#include "random_utils.h"

Partition::Partition(Graph g, int num_partitions) : graph(g), k(num_partitions), partition(graph.n, -1), nodes_in_partition(num_partitions, 0) {
    alpha = (double) graph.m * pow(k, 0.5) / pow(graph.n, gamma);
}

void Partition::stream_partition(string permutation, string status, bool enforce_balance_constraint) {
    int progress = 0;
    if (status == "progress") {
        cout << "|" << flush;
    }

    int capacity = v * (graph.n / ((float) k));
    vector<NodeID> node_ordering;
    get_node_ordering(node_ordering, permutation);

    for (int i = 0; i < graph.n; i++) {
        double max_gain = std::numeric_limits<float>::lowest();
        double gain;
        PartitionID best_partition = randomFrom(0, (int) k - 1);
        NodeID node_id = node_ordering[i];
        for (PartitionID p_id = 0; p_id < (unsigned int) k; p_id++) {
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




void Partition::print_partition_evaluation() {
    Partition::print_partition_evaluation_title();
    Partition::print_partition_evaluation(graph.n, graph.m, k, cutsize, 100 * lambda, rho);
}

const int n_width = 7;
const int m_width = 10;
const int k_width = 4;
const int lambda_width = 8;
const int rho_width = 7;

void Partition::print_partition_evaluation_title() {
    printElement("n", n_width);
    printElement("m", m_width);
    printElement("k", k_width);
    printElement("cutsize", n_width);
    printElement("lambda", lambda_width);
    printElement("rho", rho_width);
    cout << endl << "---------------------------------------------------" << endl;
}

void Partition::print_partition_evaluation(int n, int m, int num_partition, int cutsize, double lambda, double rho) {
    printElement(n, n_width);
    printElement(m, m_width);
    printElement(num_partition, k_width);
    printElement(cutsize, n_width);
    printElement(lambda, lambda_width);
    printElement(rho, rho_width);
    cout << endl;
}

void Partition::get_node_ordering(vector<NodeID> &node_ordering, string permutation) {
    if (permutation == "standard" || permutation == "random") {
        node_ordering.resize(graph.n);
        for (NodeID i = 0; i < (unsigned int) graph.n; ++i) {
            node_ordering[i] = i;
        }
        if (permutation == "random") {
            random_shuffle(node_ordering.begin(), node_ordering.end());
        }
    } else if (permutation == "dfs" || permutation == "bfs") {
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(0, graph.n - 1);
        int start_node = dis(gen);
        if (permutation == "dfs") {
            dfs(node_ordering, start_node);
        } else {
            bfs(node_ordering, start_node);
        }
    }
}

void Partition::dfs(vector<NodeID> &ordering, NodeID start_node_id) {
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

    NodeID node_id = 0;
    for (bool visited_node : already_visited) {
        if (!visited_node) {
            ordering.push_back(node_id);
        }
        node_id++;
    }
}

void Partition::bfs(vector<NodeID> &ordering, NodeID start_node_id) {
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

    NodeID node_id = 0;
    for (bool visited_node : visited) {
        if (!visited_node) {
            ordering.push_back(node_id);
        }
        node_id++;
    }
}

double Partition::fennel_gain(NodeID node_id, PartitionID p_id, int nodes_in_partition) {
    int neighbours_in_partition = 0;
    for (NodeID adj: graph.nodes[node_id].adjacents) {
        if (partition[adj] == p_id) {
            neighbours_in_partition++;
        }
    }

    double intra_partition_cost = alpha * gamma * sqrt(nodes_in_partition);
    return neighbours_in_partition - intra_partition_cost;
}

void Partition::compute_cut_size() {
    cutsize = 0;
    for (NodeID node_id = 0; node_id < (unsigned int) graph.n; node_id++) {
        for (NodeID neighbor_id : graph.nodes[node_id].adjacents) {
            if (partition[neighbor_id] != partition[node_id]) {
                cutsize++;
            }
        }
    }
    cutsize /= 2; // Divide by 2 to account for doubled edges
}

void Partition::eval_partition() {
    compute_cut_size();
    lambda = cutsize / graph.m;
    int max_load = 0;
    for (PartitionID i = 0; i < (unsigned int) k; i++) {
        int load = count(partition.begin(), partition.end(), i);
        if (load > max_load) {
            max_load = load;
        }
    }
    rho = max_load / (partition.size() / (double) k);
}

void Partition::log_progress(int i, int &progress) {
    int log_steps = 10;
    if ((i % (int) (1.0 / log_steps * graph.n)) == 0) {
        cout << "o" << flush;
        progress++;
    }
}
