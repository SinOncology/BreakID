//
// Created by jinlf on 7/31/17.
//

#include "util_cluster.h"

void init_cluster(cluster_struct &main_cluster, long distance_threshold,
                  vector<point> &points, int linkage_type)
{
  clock_t matrix_start, matrix_end, add_leaves_start, add_leaves_end, merge_start, merge_end;
  
  matrix_start = clock();
  double **matrix = generate_distance_matrix(points);
  matrix_end = clock();
//  cout << "the generate matrix step  costs time: " <<
//       (matrix_end - matrix_start) / (CLOCKS_PER_SEC * 60) << " minutes "
//       << (matrix_end - matrix_start) / CLOCKS_PER_SEC - ((matrix_end - matrix_start) / (CLOCKS_PER_SEC * 60)) * 60
//       << " seconds " << std::endl;
  cout << "the generate matrix step  costs time: " <<
       (matrix_end - matrix_start) / double(CLOCKS_PER_SEC) << " seconds " << std::endl;
  
  main_cluster.distance_matrix = matrix;
  main_cluster.num_points = points.size();
  main_cluster.num_nodes = 0;
  main_cluster.num_root_clusters = 0;
  add_leaves_start = clock();
  add_leaves(main_cluster, points);
  add_leaves_end = clock();
//  cout << "the add leaves step  costs time: " <<
//       (add_leaves_end - add_leaves_start) / (CLOCKS_PER_SEC * 60) << " minutes "
//       << (add_leaves_end - add_leaves_start) / CLOCKS_PER_SEC - ((add_leaves_end - add_leaves_start) / (CLOCKS_PER_SEC * 60)) * 60
//       << " seconds " << std::endl;
  
  cout << "the add leaves step  costs time: " <<
       (add_leaves_end - add_leaves_start) / double(CLOCKS_PER_SEC) << " seconds " << std::endl;
  
  merge_start = clock();
  merge_clusters(main_cluster, distance_threshold, linkage_type);
  merge_end = clock();
//  cout << "the merge cluster step  costs time: " <<
//       (merge_start - merge_end) / (CLOCKS_PER_SEC * 60) << " minutes "
//       << (merge_start - merge_end) / CLOCKS_PER_SEC - ((merge_start - merge_end) / (CLOCKS_PER_SEC * 60)) * 60
//       << " seconds " << std::endl;
  cout << "the merge cluster step  costs time: " <<
       (merge_end - merge_start) / double(CLOCKS_PER_SEC) << " seconds " << std::endl;
  
}

double **generate_distance_matrix(vector<point> &points)
{
  double **matrix = (double **) (calloc(points.size(), sizeof(double *)));
  if (matrix)
  {
    for (int i = 0; i < points.size(); ++i)
    {
      matrix[i] = (double *) calloc(points.size(), sizeof(double));
      if (!matrix[i])
      {
        for (i = 0; i < points.size(); ++i)
        {
          free(matrix[i]);
        }
        free(matrix);
        matrix = NULL;
        break;
      }
    }
    for (int i = 0; i < points.size(); ++i)
    {
      for (int j = 0; j < points.size(); ++j)
      {
        *(*(matrix + i) + j) = euclidean_distance(points[i].pos, points[j].pos);
      }
    }
  }
  return matrix;
}

double euclidean_distance(coordinate &pos_1, coordinate &pos_2)
{
  double dist;
  dist = sqrt(pow(pos_1.x - pos_2.x, 2) + pow(pos_1.y - pos_2.y, 2));
  return dist;
}

void add_leaves(cluster_struct &main_cluster, vector<point> &points)
{
  for (int i = 0; i < main_cluster.num_points; ++i)
  {
    add_leaf(main_cluster, points[i]);
    update_neighbours(main_cluster, 1);
  }
  
}

void add_leaf(cluster_struct &main_cluster, const point &single_point)
{
  node temp_node;
  temp_node.centroid = single_point.pos;
  temp_node.type = LEAF_NODE;
  temp_node.is_root = 1;
  temp_node.height = 0;
  temp_node.num_points = 1;
  temp_node.points.push_back(main_cluster.num_nodes);
  temp_node.label = single_point.label;
  temp_node.neighbours = NULL;
  main_cluster.nodes.push_back(temp_node);
  main_cluster.num_root_clusters++;
  main_cluster.num_nodes++;
}

void update_neighbours(cluster_struct &main_cluster, int linkage_type)
{
//  node *current_node = &(main_cluster->nodes[index]);
  int current_node_index = (int) main_cluster.nodes.size() - 1;
  if (main_cluster.nodes[current_node_index].type != NOT_USED)
  {
    int root_cluster_seen = 1;
    int target = current_node_index;
    while (root_cluster_seen < main_cluster.num_root_clusters)
    {
      target--;
      if (main_cluster.nodes[target].type == NOT_USED)
      {
        break;
      }
      if (main_cluster.nodes[target].is_root)//只add root 节点到neighbour
      {
        root_cluster_seen++;
        add_neighbour(main_cluster, current_node_index, target, linkage_type);
      }
    }
  }
}

void add_neighbour(cluster_struct &main_cluster,
                   int current_node_index, int target, int linkage_type)
{
  neighbour *current_neighbour = (neighbour *) calloc(1, sizeof(neighbour));
  if (current_neighbour)
  {
    current_neighbour->target = target;
    current_neighbour->distance = get_distance(main_cluster, current_node_index, target, linkage_type);

//    std::cerr << "single linkage distance of node No."<< index <<" and node No." << target <<" is: "<<current_neighbour->distance<<"\n";
    
    if (main_cluster.nodes[current_node_index].neighbours)
    {
      insert_sorted(current_neighbour, main_cluster, current_node_index);
    }
    else
    {
      main_cluster.nodes[current_node_index].neighbours = current_neighbour;
    }
  }
}

double get_distance(cluster_struct &main_cluster,
                    int current_node_index, int target, int linkage_type)
{
  if (current_node_index < main_cluster.num_points && target < main_cluster.num_points)
  {
//    return main_cluster.distance_matrix[current_node_index][target];
    return *(*(main_cluster.distance_matrix + current_node_index) + target);
  }
  else
  {
    node a = main_cluster.nodes[current_node_index];
    node b = main_cluster.nodes[target];
    double dist;
    switch (linkage_type)
    {
    case 1:
      dist =
        average_linkage(main_cluster.distance_matrix,
                        a.points, b.points, a.num_points, b.num_points);
      break;
    case 2:
      dist =
        complete_linkage(main_cluster.distance_matrix,
                         a.points, b.points, a.num_points, b.num_points);
      break;
    
    case 3:
      dist =
        single_linkage(main_cluster.distance_matrix,
                       a.points, b.points, a.num_points, b.num_points);
      
      
      break;
    
    default:
      dist =
        average_linkage(main_cluster.distance_matrix,
                        a.points, b.points, a.num_points, b.num_points);
    }
    return dist;
  }
}

double average_linkage(double **matrix, vector<int> a, vector<int> b,
                       int m, int n)
{
  double total = 0.0;
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
//      total += matrix[a[i]][b[j]];
      total += *(*(matrix + a[i]) + b[j]);
    }
  }
  double dist = total / (m * n);
  return dist;
}

double complete_linkage(double **matrix, vector<int> a, vector<int> b,
                        int m, int n)
{
  double d, max_dist = 0.0;
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      d = *(*(matrix + a[i]) + b[j]);
      if (d > max_dist)
        max_dist = d;
    }
  }
  return max_dist;
}

double single_linkage(double **matrix, vector<int> a, vector<int> b,
                      int m, int n)
{
  double d, min_dist = DBL_MAX;
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      d = *(*(matrix + a[i]) + b[j]);
      if (d < min_dist)
        min_dist = d;
    }
  }
  return min_dist;
}

void insert_sorted(neighbour *current_neighbour,
                   cluster_struct &main_cluster, int current_node_index)
{
  
  
  neighbour *temp;
  temp = main_cluster.nodes[current_node_index].neighbours;
  while (temp->next)
  {
    if (temp->distance >= current_neighbour->distance)
    {
      insert_before(temp, current_neighbour, main_cluster, current_node_index);
      return;
    }
    temp = temp->next;
  }
  
  if (temp->distance > current_neighbour->distance)
  {
    insert_before(temp, current_neighbour, main_cluster, current_node_index);
  }
  else
  {
    insert_after(temp, current_neighbour);
  }
  
}

void insert_before(neighbour *temp, neighbour *current_neighbour,
                   cluster_struct &main_cluster, int current_node_index)
{
  current_neighbour->next = temp;
  if (temp->prev)
  {
    temp->prev->next = current_neighbour;
    current_neighbour->prev = temp->prev;
  }
  else
  {
    main_cluster.nodes[current_node_index].neighbours = current_neighbour;
  }
  temp->prev = current_neighbour;
}

void insert_after(neighbour *temp, neighbour *current_neighbour)
{
  current_neighbour->prev = temp;
  temp->next = current_neighbour;
}

void merge_clusters(cluster_struct &main_cluster, long distance_threshold, int linkage_type)
{
  int first, second;
  double best_dist;
  while (main_cluster.num_root_clusters > 1)
  {
    best_dist = DBL_MAX;
    first = -1, second = 0;
    int return_code = find_cluster_to_merge(main_cluster, first, second, best_dist);
    if (return_code != -1 && best_dist <= distance_threshold)
    {
      merge(main_cluster, first, second, linkage_type);
    }
    else
    {
      break;
    }
  }
  
}

int &find_cluster_to_merge(cluster_struct &main_cluster, int &first, int &second, double &best_dist)
{
  int root_cluster_seen = 0;
  int j = main_cluster.num_nodes;
  while (root_cluster_seen < main_cluster.num_root_clusters)
  {//遍历所有根cluster节点
    j--;
    //非根结点就不看
    if (main_cluster.nodes[j].type == NOT_USED || !main_cluster.nodes[j].is_root) {continue;}
    root_cluster_seen++;
    find_best_distance_neighbour(main_cluster, j,
                                 best_dist, first, second);
  }
  return first;
}


void find_best_distance_neighbour(cluster_struct &main_cluster,
                                  int j, double &best_dist, int &first, int &second)
{
  neighbour *neighbours = main_cluster.nodes[j].neighbours;
  while (neighbours)
  {
    if (main_cluster.nodes[neighbours->target].is_root)
    {
      if (first == -1 || neighbours->distance < best_dist)
      {
        first = j;
        second = neighbours->target;
        best_dist = neighbours->distance;
      }
      break;
    }
    neighbours = neighbours->next;
  }
}

void merge(cluster_struct &main_cluster, int &first, int &second, int linkage_type)
{
  node temp_node;
  temp_node.merged.push_back(first);
  temp_node.merged.push_back(second);
  temp_node.num_points = main_cluster.nodes[first].num_points
                         + main_cluster.nodes[second].num_points;
  temp_node.type = MERGER;
  temp_node.is_root = 1;
  temp_node.height = -1;
  
  coordinate centroid = {.x = 0.0, .y = 0.0};
  std::string label = "node";
  int to_merge[2] = {first, second};
  for (int i = 0; i < 2; ++i)
  {
    main_cluster.nodes[*(to_merge + i)].is_root = 0;
    if (main_cluster.nodes[*(to_merge + i)].height == -1 ||
        temp_node.height < main_cluster.nodes[*(to_merge + i)].height)
    {
      temp_node.height = main_cluster.nodes[*(to_merge + i)].height;
    }
    for (int j = 0; j < main_cluster.nodes[*(to_merge + i)].num_points; ++j)
    {
      temp_node.points.push_back(main_cluster.nodes[*(to_merge + i)].points[j]);
    }
    centroid.x += main_cluster.nodes[*(to_merge + i)].num_points * main_cluster.nodes[*(to_merge + i)].centroid.x;
    centroid.y += main_cluster.nodes[*(to_merge + i)].num_points * main_cluster.nodes[*(to_merge + i)].centroid.y;
    label += ":" + main_cluster.nodes[*(to_merge + i)].label;
  }
  temp_node.centroid.x = centroid.x / temp_node.num_points;
  temp_node.centroid.y = centroid.y / temp_node.num_points;
  temp_node.height++;
  temp_node.label = label;
  temp_node.neighbours = NULL;
  main_cluster.nodes.push_back(temp_node);
  main_cluster.num_nodes++;
  main_cluster.num_root_clusters--;
  update_neighbours(main_cluster, linkage_type);
}

int print_root_nodes(cluster_struct &main_cluster)
{
//  node* temp_node = main_cluster->nodes;
  int k = 0;
  for (int i = 0; i < main_cluster.num_nodes; ++i)
  {
    
    if (main_cluster.nodes[i].is_root)
    {
      k++;
//      std::cerr << "the label of root cluster is:"
//                << temp_node->label << std::endl;
//      std::cerr << "the "<<i<<"th node, it is also the  "<<k<<"th root node:\n=====================" << std::endl;
//      print_current_node(main_cluster.nodes[i]);
    }
    else {continue;}
  }
  
  return k;
}

void print_current_node(node &temp_node)
{
  if (temp_node.num_points >= 2)
  {
    std::cerr << "the nearest distance nodes of the current root node(contains more than 2 points) is:\t" <<
              temp_node.neighbours->distance << std::endl;
//    std::cerr << "it contains " << temp_node->num_points << " points" << std::endl;
//    std::string print_str = "list as below:\n";
//    for (int i = 0; i < temp_node->num_points; ++i)
//    {
//      print_str += "the " + std::to_string(temp_node->points[i]) + "th point\n";
//    }
//    std::cerr << print_str;
  }
}