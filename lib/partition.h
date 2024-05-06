#ifndef PARTITION_H
#define PARTITION_H

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
#include "graph.h"



using namespace std;

typedef unsigned int PartitionID;


template<typename T> void printElement(T t, const int& width);



class Partition {
public:
    Graph graph;
    int k;
    vector<PartitionID> partition;
    vector<int> nodes_in_partition;

    double v = 1.1;
    double gamma = 1.5;
    double alpha, lambda, rho;
    double cutsize;

    Partition(Graph g, int num_partitions);

    void stream_partition(string permutation, string status="progress", bool enforce_balance_constraint=true);

    void buffered_stream_partition(string permutation, int buffer_size, string status="progress", bool enforce_balance_constraint=true);

    void print_partition_evaluation();

    static void print_partition_evaluation_title();

    static void print_partition_evaluation(int n, int m, int num_partition, int cutsize, double lambda, double rho);

private:
    void get_node_ordering(vector<NodeID>& node_ordering, string permutation);

    void dfs(vector<NodeID>& ordering, NodeID start_node_id);

    void bfs(vector<NodeID>& ordering, NodeID start_node_id);

    double fennel_gain(NodeID node_id, PartitionID p_id, int nodes_in_partition);

    void compute_cut_size();

    void eval_partition();

    void log_progress(int i, int& progress);
};

template<typename T> void printElement(T t, const int& width) {
    const char separator = ' ';
    cout << left << setw(width) << setfill(separator) << setprecision(3) << t << "|";
}

#endif // PARTITION_H
