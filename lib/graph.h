#ifndef GRAPH_HEADER_H
#define GRAPH_HEADER_H

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
#include <fstream>
#include <sstream>

using namespace std;

typedef unsigned int NodeID;
typedef unsigned int EdgeID;

struct Node {
    vector<NodeID> adjacents;
};

class Graph {
public:
    int n, m;
    vector<Node> nodes;
    double gamma, alpha, v;

    Graph() {}

    int read_edge_list_from_file(string filename);
    int read_graph_from_metis_file(string filename);
    static Graph generate_hp_random_graph(int num_nodes, int num_partitions, double intra_prob, double inter_prob);
};

#endif // GRAPH_HEADER_H
