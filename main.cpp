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


bool ORDERING_STANDARD = true;
bool ORDERING_RANDOM = true;
bool ORDERING_DFS = true;
bool ORDERING_BFS = true;

int main(int argc, char* argv[]) {


    if (argc != 3) {
            cerr << "Usage: " << argv[0] << " <filename>" << "(optional) <k>" << endl;
            return 1;
    } else {
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
