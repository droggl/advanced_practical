#include "partition.h"
#include "../global.h"
#include "random_utils.h"
#include <list>

long getMaxRSS() {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // The maximum resident set size is in kilobytes
        return usage.ru_maxrss;
    } else {
        std::cerr << "Error getting resource usage information." << std::endl;
        // Return a sentinel value or handle the error in an appropriate way
        return -1;
    }
}

Partition::Partition(Graph g, int num_partitions) : graph(g), k(num_partitions), partition(graph.n, -1), nodes_in_partition(num_partitions, 0) {
    alpha = (double)graph.m * pow(k, 0.5) / pow(graph.n, gamma);
    capacity = v * (graph.n / ((float)k));
}

void Partition::stream_partition(vector<NodeID> &node_ordering, bool delay_nodes_enabled, unsigned eval_delay_after, string status) {
    auto start = std::chrono::high_resolution_clock::now();

    int progress = 0;
    if (status == "progress") {
        cout << "|" << flush;
    }

    list<NodeID> delayed_nodes;
    int num_delayed_nodes = 0;
    vector<int> should_be_delayed_counts(graph.n, 0);

    unsigned eval_delayed_nodes_counter = 0;
    for (unsigned i = 0; i < graph.n; i++) {
        NodeID node_id = node_ordering[i];

        if (delay_nodes_enabled) {
            if (eval_delayed_nodes_counter >= eval_delay_after - 1) { // graph.n / 10 )  { // >=  graph.n / 1000) {
                eval_delayed_nodes_counter = 0;
                for (auto it = delayed_nodes.begin(); it != delayed_nodes.end();) {
                    NodeID delayed_node_id = *it;
                    should_be_delayed_counts[delayed_node_id]++;
                    if (!should_be_delayed(delayed_node_id)) {
                        partition_node(delayed_node_id);
                        it = delayed_nodes.erase(it); // Remove node and set iterator to next element
                    } else {
                        ++it;
                    }
                }
            }
            eval_delayed_nodes_counter++;

            // Check if assignment of current node should be delayed
            if (i != 0 && should_be_delayed(node_id)) {
                num_delayed_nodes++;
                delayed_nodes.push_back(node_id);
                continue;
            }
        }

        partition_node(node_id);

        if (status == "progress") {
            log_progress(i, progress);
        }
    }

    if (delay_nodes_enabled) {
        for (NodeID delayed_node_id : delayed_nodes) {
            partition_node(delayed_node_id);
        }

        // cout << "number of delayed nodes: " << num_delayed_nodes << endl;
    }

    if (status == "progress") {
        cout << "|\n" << endl;
    }

    // for (int cnt : should_be_delayed_counts) {
    //     cout << cnt << endl;
    // }

    auto duration = chrono::high_resolution_clock::now() - start;
    t_partition_ms = chrono::duration_cast<chrono::milliseconds>(duration).count();

    maxRSS = getMaxRSS();
    eval_partition();
}

bool Partition::should_be_delayed(NodeID node_id) {
    // Do not buffer nodes with degree > MAX_DEGREE_DELAY
    if (graph.nodes[node_id].adjacents.size() > MAX_DEGREE_DELAY) {
        return false;
    }

    float percentage_of_neighbours_partitioned = 0;
    for (NodeID adj_id : graph.nodes[node_id].adjacents) {
        if (partition[adj_id] != -1) {
            percentage_of_neighbours_partitioned++;
        }
    }
    percentage_of_neighbours_partitioned /= graph.nodes[node_id].adjacents.size();
    return percentage_of_neighbours_partitioned < SHOULD_BE_DELAYED_BARRIER;
}

void Partition::partition_node(NodeID node_id) {
    double max_gain = std::numeric_limits<float>::lowest();
    double gain;
    PartitionID best_partition = randomFrom(0, (int)k - 1);
    for (PartitionID p_id = 0; p_id < k; p_id++) {
        // Check if new node in partition would exceed threshold v
        if (nodes_in_partition[p_id] >= capacity) {
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
}
//
void Partition::buffered_stream_partition(vector<NodeID> &node_ordering, int buffer_size, string status, bool delay_nodes_enabled) {
    auto start = std::chrono::high_resolution_clock::now();

    list<NodeID> delayed_nodes;
    int num_delayed_nodes = 0;

    for (NodeID start_node_id = 0; start_node_id < graph.n; start_node_id += buffer_size) {
        if (delay_nodes_enabled) {
            for (auto it = delayed_nodes.begin(); it != delayed_nodes.end();) {
                NodeID delayed_node_id = *it;
                if (!should_be_delayed(delayed_node_id)) {
                    partition_node(delayed_node_id);
                    it = delayed_nodes.erase(it); // Remove node and set iterator to next element
                } else {
                    ++it;
                }
            }
        }

        vector<NodeID> local_to_global;
        Graph subgraph = graph.extract_subgraph(node_ordering, start_node_id, buffer_size, local_to_global);
        // cout << "SUBGRAPH: " << subgraph.n << ",  " << subgraph.nodes.size() << endl;

        vector<NodeID> subgraph_bfs_ordering;
        // Partition::bfs(subgraph, subgraph_bfs_ordering, 0);
        average_bfs_depth += get_node_ordering(subgraph, subgraph_bfs_ordering, ORDERING_BFS);

        for (NodeID local_node_id : subgraph_bfs_ordering) {
            NodeID node_id = local_to_global[local_node_id];

            if (delay_nodes_enabled) {
                if (start_node_id != 0 && should_be_delayed(node_id)) {
                    num_delayed_nodes++;
                    delayed_nodes.push_back(node_id);
                    continue;
                }
            }

            partition_node(node_id);
        }
    }

    if (delay_nodes_enabled) {
        for (NodeID delayed_node_id : delayed_nodes) {
            partition_node(delayed_node_id);
        }
    }

    auto duration = chrono::high_resolution_clock::now() - start;
    t_partition_ms = chrono::duration_cast<chrono::milliseconds>(duration).count();

    if (SHOW_MAX_RSS)
        maxRSS = getMaxRSS();

    if (SHOW_AVG_BFS_DEPTH) {
        float avg_bfs_depth = average_bfs_depth / (ceil(graph.n / buffer_size));
        cout << left << setw(15) << setfill(' ') << setprecision(6) << avg_bfs_depth << "|";
    }

    eval_partition();

}

// Priority is the ratio of partitioned neighbours to total neighbours
float Partition::get_priority(NodeID node_id) {
    float number_of_partitioned_neighnours = 0;
    for (NodeID adj_id : graph.nodes[node_id].adjacents) {
        if (partition[adj_id] != -1) {
            number_of_partitioned_neighnours++;
        }
    }
    return number_of_partitioned_neighnours / graph.nodes[node_id].adjacents.size();
}

// Update the priority value of the neighbours of the node that was just partitioned in the priority queue (insert with new value)
void Partition::update_neighbours_priority(NodeID node_id, priority_queue<pq_node, vector<pq_node>, ComparePriority> &pq, vector<float> &node_id_to_priority) {
    for (NodeID adj_id : graph.nodes[node_id].adjacents) {
        // If neightbour is not partitioned and already in the priority queue, update the priority value
        if (partition[adj_id] == -1 && node_id_to_priority[adj_id] != -1) {
            // float new_priority = get_priority(adj_id);
            float new_priority = node_id_to_priority[adj_id] + 1 / (float) graph.nodes[adj_id].adjacents.size();
            node_id_to_priority[adj_id] = new_priority;
            if (new_priority == 1) {
                // If all neighbours are partitioned, remove from priority queue
                partition_node(adj_id);
                update_neighbours_priority(adj_id, pq, node_id_to_priority);
            } else {
                pq.push(pq_node(adj_id, new_priority));
            }
        }
    }
}

// Stream partitioning with priority queue
void Partition::stream_partition_with_pq(vector<NodeID> &node_ordering) {

    auto start = std::chrono::high_resolution_clock::now();

    priority_queue<pq_node, vector<pq_node>, ComparePriority> pq_delayed_nodes;
    vector<float> node_id_to_priority(graph.n, -1); // if -1: not in priority queue, else: index in priority queue with corresponding priority

    // First, assign nodes to the priority queue with priority based on the already partitioned neighbors
    for (NodeID node_id : node_ordering) {
        // If degree of node is bigger than MAX_DEGREE_DELAY, partition node directly
        bool should_be_partitioned_directly = graph.nodes[node_id].adjacents.size() > MAX_DEGREE_DELAY;
        if (should_be_partitioned_directly) {
            // Partition node and update priority of neighbours in pq
            partition_node(node_id);
            update_neighbours_priority(node_id, pq_delayed_nodes, node_id_to_priority);
        } else if (node_id_to_priority[node_id] == -1) { // If node is not already in priority queue
            // Add node to priority queue
            float priority = get_priority(node_id);
            pq_delayed_nodes.push(pq_node(node_id, priority));
            node_id_to_priority[node_id] = priority;
        }
    }
    int cnt = 0;
    // Second, empty the priority queue and assign nodes to partitions
    while (!pq_delayed_nodes.empty()) {
        pq_node pq_node = pq_delayed_nodes.top();
        pq_delayed_nodes.pop();
        if (pq_node.priority == node_id_to_priority[pq_node.id]) { // If priority is the current priority of the node
            // Partition node and update priority of neighbours in pq
            partition_node(pq_node.id);
            update_neighbours_priority(pq_node.id, pq_delayed_nodes, node_id_to_priority);
        } else {
            cnt++;
        }
    }

    auto duration = chrono::high_resolution_clock::now() - start;
    t_partition_ms = chrono::duration_cast<chrono::milliseconds>(duration).count();
    maxRSS = getMaxRSS();
    eval_partition();
}

const int ordering_width = 10;
const int configuration_width = 25;
const int n_width = 7;
const int m_width = 10;
const int k_width = 4;
const int lambda_width = 8;
const int rho_width = 7;
const int t_width = 6;
const int maxRSS_width = 12;

void Partition::print_partition_evaluation_title() {
    printElement("ordering", ordering_width);
    printElement("configuration", configuration_width);
    // printElement("n", n_width);
    // printElement("m", m_width);
    // printElement("k", k_width);
    printElement("cutsize", n_width);
    printElement("lambda", lambda_width);
    printElement("rho", rho_width);
    printElement("t (ms)", t_width);
    if (SHOW_MAX_RSS)
        printElement("maxRSS (KB)", maxRSS_width);
    cout << endl
         << "-----------------------------------------------------------------------" << endl;
}

void Partition::print_partition_evaluation(string configuration, string ordering, bool print_title) {
    if (print_title) {
        Partition::print_partition_evaluation_title();
    }
    Partition::print_partition_evaluation(configuration, ordering, graph.n, graph.m, k, cutsize, 100 * lambda, rho, t_partition_ms, maxRSS);
}

void Partition::print_partition_evaluation(string configuration, string ordering, int n, int m, int num_partition, int cutsize, double lambda, double rho, int64_t t_partition_ms, long maxRSS) {
    printElement(ordering, ordering_width);
    printElement(configuration, configuration_width);
    // printElement(n, n_width);
    // printElement(m, m_width);
    // printElement(num_partition, k_width);
    printElement(cutsize, n_width);
    printElement(lambda, lambda_width);
    printElement(rho, rho_width);
    printElement(t_partition_ms, t_width);
    if (SHOW_MAX_RSS)
        printElement(maxRSS, maxRSS_width);
    cout << endl;
}

float Partition::get_node_ordering(Graph &g, vector<NodeID> &node_ordering, string permutation) {
    float avg_bfs_depth = 0;
    if (permutation == ORDERING_NATURAL || permutation == ORDERING_RANDOM) {
        node_ordering.resize(g.n);
        for (NodeID i = 0; i < g.n; ++i) {
            node_ordering[i] = i;
        }
        if (permutation == ORDERING_RANDOM) {
            random_shuffle(node_ordering.begin(), node_ordering.end());
        }
    } else if (permutation == ORDERING_DFS || permutation == ORDERING_BFS) {
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(0, g.n - 1);
        int start_node = dis(gen);
        if (permutation == ORDERING_DFS) {
            dfs(g, node_ordering, start_node);
        } else {
            vector<bool> visited(g.nodes.size(), false);
            bfs(g, node_ordering, start_node, visited);

            NodeID node_id = 0;
            int cnt = 1;
            // cout << "checking unvisited nodes, visited size: " << visited.size() << " number: " << count(visited.begin(), visited.end(), false) << endl;
            for (bool visited_node : visited) {
                // cout << "node_id: " << node_id << " " << visited_node << endl;
                if (!visited_node) {
                    cnt++;
                    // bfs(g, node_ordering, node_id, visited);
                    node_ordering.push_back(node_id);
                }
                node_id++;
            }
            avg_bfs_depth += g.n / (float)cnt;
            // cout << "average bfs depth: " << g.n / (float) cnt << ", number of BFSs: " << cnt << endl;
            // cout << "NUMBER OF BFSs PERFORMED: " << cnt << endl;
            // cout << "n: " << g.n << ", node ordering size: " << node_ordering.size() << endl;
        }
    }
    return avg_bfs_depth;
}

void Partition::dfs(Graph &g, vector<NodeID> &ordering, NodeID start_node_id) {
    vector<bool> already_visited(g.n, false);
    stack<NodeID> stack;
    stack.push(start_node_id);

    while (!stack.empty()) {
        NodeID current_node_id = stack.top();
        stack.pop();

        if (!already_visited[current_node_id]) {
            already_visited[current_node_id] = true;
            ordering.push_back(current_node_id);

            for (NodeID adjacent_id : g.nodes[current_node_id].adjacents) {
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

void Partition::bfs(Graph &g, vector<NodeID> &ordering, NodeID start_node_id, vector<bool> &visited) {
    queue<NodeID> queue;

    queue.push(start_node_id);
    visited[start_node_id] = true;
    while (!queue.empty()) {
        NodeID current_node_id = queue.front();
        queue.pop();

        ordering.push_back(current_node_id);
        for (NodeID adjacent_id : g.nodes[current_node_id].adjacents) {
            visited[adjacent_id];
            if (!visited[adjacent_id]) {
                queue.push(adjacent_id);
                visited[adjacent_id] = true;
            }
        }
    }
}

double Partition::fennel_gain(NodeID node_id, PartitionID p_id, int nodes_in_partition) {
    int neighbours_in_partition = 0;
    for (NodeID adj : graph.nodes[node_id].adjacents) {
        if (partition[adj] == p_id) {
            neighbours_in_partition++;
        }
    }

    double intra_partition_cost = alpha * gamma * sqrt(nodes_in_partition);
    return neighbours_in_partition - intra_partition_cost;
}

void Partition::compute_cut_size() {
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

void Partition::eval_partition() {
    compute_cut_size();
    lambda = cutsize / graph.m;
    int max_load = 0;
    for (PartitionID i = 0; i < k; i++) {
        int load = count(partition.begin(), partition.end(), i);
        if (load > max_load) {
            max_load = load;
        }
    }
    rho = max_load / (partition.size() / (double)k);
}

void Partition::log_progress(int i, int &progress) {
    int log_steps = 10;
    if ((i % (int)(1.0 / log_steps * graph.n)) == 0) {
        cout << "o" << flush;
        progress++;
    }
}
