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
#include "lib/graph.h"
#include "lib/partition.h"

using namespace std;

unordered_map<string, bool> ORDERING_ENABLED = {
    {ORDERING_DFS,      false},
    {ORDERING_BFS,      true},
    {ORDERING_RANDOM,   false},
    {ORDERING_NATURAL,  true}
};

int main(int argc, char* argv[]) {


    if (argc != 4) {
            cerr << "Usage: " << argv[0] << " <filename>" << "(optional) <k>" << endl;
            return 1;
    } else {
        int k =  atoi(argv[2]);
        unsigned batch_size = atoi(argv[3]);

        string filename = argv[1];


        Graph graph = Graph();
        int ret_val = graph.read_graph_from_metis_file(filename);
        if (ret_val != 0) {
            return ret_val;
        }

        cout << endl << "# ---- Graph info ----- #" << endl;
        cout << "n: " << graph.n << ", m:" << graph.m << ", k: " << k << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        vector<NodeID> node_ordering;
        Partition p = Partition(graph, k);
        for (const auto& pair : ORDERING_ENABLED) {
            const string& ordering = pair.first;
            bool enabled = pair.second;

            if (enabled) {
                cout << endl << "Ordering: " << ordering << endl;
                Partition::print_partition_evaluation_title();

                node_ordering.clear();
                Partition::get_node_ordering(graph, node_ordering, ordering);

                p = Partition(graph, k);
                p.stream_partition(node_ordering, false);
                p.print_partition_evaluation("streaming");

                p = Partition(graph, k);
                p.stream_partition(node_ordering, true, 1);
                p.print_partition_evaluation("streaming+delay(1)");

                p = Partition(graph, k);
                p.stream_partition(node_ordering, true, graph.n/100);
                p.print_partition_evaluation("streaming+delay(n/100)");

                p = Partition(graph, k);
                p.buffered_stream_partition(node_ordering, batch_size, "none", false);
                p.print_partition_evaluation("subgraph" + to_string(batch_size));

                p = Partition(graph, k);
                p.buffered_stream_partition(node_ordering, batch_size, "none", true);
                p.print_partition_evaluation("subgraph"+to_string(batch_size)+"+delay");
                cout << "_________________________________________________________" << endl;
            }
        }


        cout << endl;

        auto duration = chrono::high_resolution_clock::now() - begin;
        auto ms = chrono::duration_cast<chrono::milliseconds>(duration).count();
        cout << endl << "Running time: " << ms << "ms" << endl << endl;
    }


    return 0;
}
