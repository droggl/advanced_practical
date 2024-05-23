#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <unordered_map>
#include "lib/graph.h"
#include "lib/partition.h"
#include "global.h"


unordered_map<string, bool> ORDERING_ENABLED = {
    {ORDERING_DFS,      false},
    {ORDERING_BFS,      false},
    {ORDERING_RANDOM,   true},
    {ORDERING_NATURAL,  true}
};

map<string, vector<NodeID>> saved_node_orderings;

void get_cached_node_ordering(Graph& g, vector<NodeID>& node_ordering, string ordering) {
    if (saved_node_orderings.find(ordering) == saved_node_orderings.end()) {
        node_ordering.clear();
        Partition::get_node_ordering(g, node_ordering, ordering);
        saved_node_orderings[ordering] = node_ordering;
    } else {
        node_ordering = saved_node_orderings[ordering];
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
            cerr << "Usage: " << argv[0] << " <filename>" << "(optional) <k>" << endl;
            return 1;
    } else {
        int k =  atoi(argv[2]);
        unsigned batch_size = atoi(argv[3]);

        string filename = argv[1];

        Graph graph = Graph();
        auto begin_reading_graph = std::chrono::high_resolution_clock::now();
        int ret_val = graph.read_graph_from_metis_file(filename);
        if (ret_val != 0) {
            return ret_val;
        }
        auto duration_reading_graph = chrono::high_resolution_clock::now() - begin_reading_graph;
        auto ms_reading_graph = chrono::duration_cast<chrono::milliseconds>(duration_reading_graph).count();
        cout << endl << "Reading graph t: " << ms_reading_graph << "ms" << endl << endl;

        cout << endl << "# ---- Graph info ----- #" << endl;
        cout << "n: " << graph.n << ", m:" << graph.m << ", k: " << k << endl << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        vector<NodeID> node_ordering;
        Partition p = Partition(graph, k);

        Partition::print_partition_evaluation_title();
        for (string ordering : {ORDERING_NATURAL, ORDERING_RANDOM, ORDERING_BFS}) {
            get_cached_node_ordering(graph, node_ordering, ordering);
            p = Partition(graph, k);
            p.stream_partition(node_ordering, false);
            p.print_partition_evaluation("streaming", ordering);
        }

        vector<int> batch_sizes; // = {16384, 32768}; // {1, 10, 64, 512, 1024, 2048, 4096, 8192, 16384, 32768};

        int j = 0;
        for (int i = graph.n; i > 100; i /= 2) {
            if (j<2) {
                batch_sizes.push_back(i);
            } else if (j % 2 == 0) {
                batch_sizes.push_back(i);
            }
            j++;
        }
        sort(batch_sizes.begin(), batch_sizes.end());

        cout << endl;
        if(SHOW_AVG_BFS_DEPTH)
            cout << left << setw(15) << setfill(' ') << "AVG BFS depth" << "|";
        Partition::print_partition_evaluation_title();
        for (const auto& pair : ORDERING_ENABLED) {
            const string& ordering = pair.first;
            bool enabled = pair.second;

            if (enabled) {
                get_cached_node_ordering(graph, node_ordering, ordering);

                for(int delay_after : batch_sizes) {
                    p = Partition(graph, k);
                    p.stream_partition(node_ordering, true, delay_after);
                    if(SHOW_AVG_BFS_DEPTH)
                        cout << left << setw(15) << setfill(' ') << setprecision(6) << p.average_bfs_depth << "|";
                    p.print_partition_evaluation("streaming+delay(" + to_string(delay_after) + ")", ordering);
                }

                for (int batch_size_new : batch_sizes) {
                    p = Partition(graph, k);
                    p.buffered_stream_partition(node_ordering, batch_size_new, "none", false);
                    p.print_partition_evaluation("subgraph" + to_string(batch_size_new), ordering);
                }

                for (int batch_size_new : batch_sizes) {
                    p = Partition(graph, k);
                    p.buffered_stream_partition(node_ordering, batch_size_new, "none", true);
                    p.print_partition_evaluation("subgraph"+to_string(batch_size_new)+"+delay", ordering);
                }

                cout << "______________________________________________________________________" << endl;
            }
        }


        cout << endl;

        auto duration = chrono::high_resolution_clock::now() - begin;
        auto ms = chrono::duration_cast<chrono::milliseconds>(duration).count();
        cout << endl << "Running time: " << ms << "ms" << endl << endl;
    }


    return 0;
}
