//
// Created by jinlf on 7/31/17.
//
#pragma once
#ifndef SINOTOOLS_UTIL_CLUSTER_C_H
#define SINOTOOLS_UTIL_CLUSTER_C_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <list>

#endif //SINOTOOLS_UTIL_CLUSTER_C_H
#define NOT_USED 0
#define LEAF_NODE 1
#define MERGER 2

using namespace std;
typedef struct coordinate coordinate;
typedef struct point point;
typedef struct node node;
typedef struct neighbour neighbour;
typedef struct cluster_struct cluster_struct;


typedef struct coordinate
{
  double x, y;
} coordinate;

typedef struct point
{
  coordinate pos;
  std::string label;
  int cluster_id=-1;
} point;

typedef struct neighbour
{
  int target;
  double distance;
  neighbour *prev;
  neighbour *next;
} neighbour;

typedef struct node
{
  int type;
  int is_root;
  int height;
  coordinate centroid;
  string label;
  vector<int> merged;
  int num_points;
  vector<int> points;
  neighbour *neighbours;
} node;

typedef struct cluster_struct
{
  unsigned long num_points;
  int num_root_clusters;
  int num_nodes;
  vector<node> nodes;
  double** distance_matrix;
} cluster_struct;


void init_cluster(cluster_struct &, long, vector<point> &, int);

double ** generate_distance_matrix(vector<point> &points);

double euclidean_distance(coordinate &pos_1, coordinate &pos_2);

void add_leaves(cluster_struct &, vector<point> &);

void add_leaf(cluster_struct &main_cluster, const point &single_point);

void update_neighbours(cluster_struct &main_cluster, int linkage_type);

void add_neighbour(cluster_struct &main_cluster,
                   int current_node_index, int target, int linkage_type);

double get_distance(cluster_struct &main_cluster,
                    int current_node_index, int target, int linkage_type);

double average_linkage(double **matrix, vector<int> a, vector<int> b,
                       int m, int n);

double complete_linkage(double **matrix, vector<int> a, vector<int> b,
                        int m, int n);

double single_linkage(double **matrix, vector<int> a, vector<int> b,
                      int m, int n);

void insert_sorted(neighbour *current_neighbour,
                   cluster_struct &main_cluster, int current_node_index);

void insert_before(neighbour *temp, neighbour *current_neighbour,
                   cluster_struct &main_cluster, int current_node_index);

void insert_after(neighbour *temp, neighbour *current_neighbour);

void merge_clusters(cluster_struct &main_cluster, long distance_threshold, int linkage_type);

int &find_cluster_to_merge(cluster_struct &main_cluster,
                          int &first, int &second, double &best_dist);

void find_best_distance_neighbour(cluster_struct &main_cluster,
                                  int j, double &best_dist, int &first, int &second);

void merge(cluster_struct &main_cluster, int &first,int &second, int linkage_type);

int print_root_nodes(cluster_struct &main_cluster);

void print_current_node(node &temp_node);
