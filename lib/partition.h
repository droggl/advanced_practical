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

typedef int PartitionID;

enum Ordering {NATURAL, RANDOM, BFS, DFS};

const string ORDERING_NATURAL = "natural";
const string ORDERING_RANDOM = "random";
const string ORDERING_BFS = "bfs";
const string ORDERING_DFS = "dfs";

float SHOULD_BE_DELAYED_BARRIER = 0.15;


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
    int capacity;

    int64_t t_partition_ms;

    Partition(Graph g, int num_partitions);

    void stream_partition(vector<NodeID>& node_ordering, bool delayed_nodes_enabled, unsigned eval_delay_after=1, string status="none");

    void buffered_stream_partition(vector<NodeID>& node_ordering, int buffer_size, string status, bool delay_nodes_enabled=true);

    void print_partition_evaluation(string configuration, bool print_title=false);

    static void print_partition_evaluation_title();

    static void print_partition_evaluation(string configuration, int n, int m, int num_partition, int cutsize, double lambda, double rho, int64_t t_partition_ms);

    static void get_node_ordering(Graph& g, vector<NodeID>& node_ordering, string permutation);
private:
    PartitionID partition_node(NodeID node_id);

    bool should_be_delayed(NodeID node_id);

    static void dfs(Graph& g, vector<NodeID>& ordering, NodeID start_node_id);

    static void bfs(Graph& g, vector<NodeID>& ordering, NodeID start_node_id);

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
