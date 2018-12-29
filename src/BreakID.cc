
#include <algorithm>
#include <zlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <sstream>
#include <htslib/sam.h>
#include <sam.h>

#include "util_bam.h"
#include "util_bed.h"
#include "util_cluster.h"
#include "installdir.h"
#include "RefSeqTranscript.h"
#include "BamAlignment.h"
#include <getopt.h>
#include <cmath>

using namespace std;

string enspan_help = "      SYNOPSIS \n \t sinotools enspan <options>\n\n \
     DESCRIPTION\n \
     \t -h -? -help \t help\n \
     \t -i        \t input bam-file\n \
     \t -o        \t output file (prefix only)\n \
     \t -q        \t encompassing reads quality thresholds  [20]\n\
     \t -b        \t genome build (hg18,hg19) [hg19]\n \
     \t -t        \t distance relative to (sqrt(2)*(insert size mean +3* insert size sd))  [2]\n \
     \t -fast     \t use the fast cluster strategy [default no] \n \
     \t -all       \t no filter enspan out [default is filter]  \n ";

//structures


typedef struct
{
  string id;
  string qname;
  long p1_flag;
  long p2_flag;
  string p1_chr;
  string p2_chr;
  uint32_t p1_pos;
  uint32_t p2_pos;
  long p1_mapq;
  long p2_mapq;
  char p1_strand;
  char p2_strand;
  int is_isolated;
  int cluster;
  string cluster_id;
  uint32_t p1_chr_pos;
  uint32_t p2_chr_pos;
} discordant_pair;

typedef struct
{
  long id;
  string p1_chr;
  uint64_t p1_mean_pos;
  uint32_t p1_min_pos;
  uint32_t p1_max_pos;
  uint32_t p1_exact_pos;
  string p2_chr;
  uint64_t p2_mean_pos;
  uint32_t p2_min_pos;
  uint32_t p2_max_pos;
  int32_t p2_exact_pos;
  long n_split_read;
  long n_discordant_pair;
  bool inv;
  string inv_type = ".";
  string fusion_type = ".";
  string discordant_reads;
  string split_reads;
  string p1_behalf_gene;
  string p1_gene_part = "0";
  string p1_genes;
  string p2_behalf_gene;
  string p2_gene_part = "0";
  string p2_genes;
  string p1_rpt;
  string p2_rpt;
  string p1_strand;
  string p2_strand;
  string p1_exon_info;
  string p2_exon_info;
  string p1_part = "";
  string p2_part = "";
  bool hotspot = false;
  bool cosmic = false;
  set<string> drp_type_set;
//  string cosmicID = "";
//  string fusionID = "";
  string fusion_pair = "";
  string p1_bp_exon = "";
  string p2_bp_exon = "";
  string up_gene = "";
  string down_gene = "";
  bool sino_pair_match = false;
  bool cosmic_pair_match = false;
  double p1_bp_depth = 0;
  double p2_bp_depth = 0;
  double p1_coverage = 0;
  double p2_coverage = 0;
  float p1_alle_freq = 0.0;
  float p2_alle_freq = 0.0;
  bool is_rpt = false;
} cluster_info;


typedef struct
{
  string primary_chr = "";
  string secondary_chr = "";
  uint32_t primary_start = 0;
  uint32_t secondary_start = 0;
  uint32_t primary_end = 0;
  uint32_t secondary_end = 0;
  string primary_cigar_str = "";
  string secondary_cigar_str = "";
  string sa_tag = "";
  uint32_t primary_bp = 0;
  uint32_t secondary_bp = 0;
  bool secondary = false;
  int flag = -1;
  string read_name = "";
  bam1_t *current_align;
} split_align_pair;

typedef struct
{
  int32_t p1_pos = -1;
  int32_t p2_pos = -1;
  string p1_match_part = "";
  string p2_match_part = "";
} bp_pos_pair;


typedef struct
{
  string p1_chr = "";
  string p2_chr = "";
  int32_t p1_bp = -1;
  int32_t p2_bp = -1;
  string p1_part = "";
  string p2_part = "";
  int encompass_num = 0;
} breakpoint_pair;


//prototypes
void
scan_discordant_pairs(const string &inp_file, const string &build, long qual, double w,
                      map<string, vector<discordant_pair>> &enspan_map, string nib_dir);

void annotate_cluster_for_sa_tag(vector<cluster_info> &clusters, string nib_dir);

void get_mean_insert_size(string input_bam, vector<double> &insert);

void get_chr_pair_set(vector<discordant_pair> &enspan, set<string> &chr_pair_set);

void mask_pairs_chr_pos(vector<discordant_pair> &enspan, long distance);

void build_pair_array(vector<discordant_pair> &, vector<point> &);

bool cmp_p1_enspan_pairs(discordant_pair a, discordant_pair b)
{
  return (a.p1_chr_pos < b.p1_chr_pos);
}

bool cmp_p2_enspan_pairs(discordant_pair a, discordant_pair b)
{
  return (a.p2_chr_pos < b.p2_chr_pos);
}

bool cmp_enspan_id(discordant_pair a, discordant_pair b)
{
  return (a.cluster < b.cluster);
}

bool cmp_cluster(cluster_info a, cluster_info b)
{
  return (a.n_discordant_pair > b.n_discordant_pair);
}

int
find_cluster_pairs_enspan_ahc(vector<discordant_pair> &enspan, double distance_threshold, int distance_type,
                              int min_reads_per_cluster);

void
add_cluster_id_for_enspan_vec(cluster_struct &main_cluster, vector<discordant_pair> &enspan, int min_reads_per_cluster);

void add_enspan_point_id(vector<discordant_pair> &enspan_vec);

void remove_isolated_pairs(vector<discordant_pair> &enspans, double w);

void write_enspan_out(string out_file, vector<cluster_info> &cluster, bool filter);

void write_enspan_params(string inp_file, string out_file, string build, double w, long qual);

void
add_exon_anno(vector<RefSeqTranscript> &txpts, string chr1, string chr2, long &p1_pos, long &p2_pos,
              vector<string> &exon_infos, const string &p1_part, const string &p2_part);

void add_exon_num_anno(RefSeqTranscript &txpt, long &pos, string &chr, vector<int> &exon_nums);

int find_cluster_pairs_enspan_fast(vector<discordant_pair> &enspan, double w, int min_reads);

void
findClusterBreakPointInfoSaTag(string bam_file, vector<discordant_pair> &enspan, double w,
                               vector<cluster_info> &cluster_vec, vector<bam1_t *> &split_reads, string nib_dir);

void find_sa_reads(samfile_t *fp, const string region_chr, uint32_t region_start, uint32_t region_end,
                   map<string, vector<split_align_pair>> &encompassing_map, bam_index_t *idx_bam);

void findEncompassingReadsAndBreakPointInfo(map<long, cluster_info> &cluster_map, const string bam_file, const int w,
                                            vector<bam1_t *> &split_reads);

void find_bp_pair(map<string, vector<split_align_pair>> &p1_encompass_map,
                  map<string, vector<split_align_pair>> &p2_encompass_map, breakpoint_pair &bp_pair,
                  const string &p1_chr, const string &p2_chr, vector<bam1_t *> &split_reads, int bp_pos_error);

string determine_fusion_type_from_drp(cluster_info &cluster);

//main program
int main(int argc, char *argv[])
{
  clock_t start, end;
  clock_t scan_start, scan_end;
  clock_t cluster_start, cluster_end;
  clock_t bp_find_start=0, bp_find_end=0;
  start = clock();
  int longindex, opt;
  // option definition
  static struct option longopts[] = {
    {"help", 0, 0, 'h'},
    {"help", 0, 0, '?'},
    {"i",    1, 0, 1},
    {"o",    1, 0, 2},
    {"q",    1, 0, 3},
    {"n",    1, 0, 4},
    {"b",    1, 0, 5},
    {"fast", 0, 0, 6},
    {"t",    0, 0, 7},
    {"all",  0, 0, 0}
//    {0,       0, 0, 0}
  };
  string inp_file = "";
  string out_file = "";
  int qual = 20;
  int times = 2;
  string build = "hg19";
  
  int distance_type = 1; //add by jin
  int min_reads_per_cluster = 2;
  bool fast_cluster = false;
 
  optind = 0;
  bool filter = true;
  string nib_dir = "";
  //parse command line arguments
  while ((opt = getopt_long_only(argc, argv, "h?", longopts, &longindex)) != -1)
  {
    switch (opt)
    {
    case 'h':
      cerr << enspan_help;
      exit(1);
    case '?':
      cerr << enspan_help;
      exit(1);
    case 1:
      inp_file = (string) optarg;
      break;
    case 2:
      out_file = (string) optarg;
      break;
    case 3:
      qual = (int) abs(atol(optarg));
      break;
    case 4:
      nib_dir = (string) optarg;
      break;
    case 5:
      build = (string) optarg;
      if (build != "hg18" && build != "hg19")
      {
        cerr << "Error: please use hg18, hg19. Falling back to default value: hg19.\n";
        build = "hg19";
      }
      break;
    case 6:
      fast_cluster = true;
      break;
    case 7:
      times = (int) abs(atol(optarg));
      break;
    case 0:
      filter = false;
      break;
    default:
      cerr << "Error: cannot parse arguments.\n";
      exit(1);
    }
  }
  
  if (inp_file.empty() || out_file.empty())
  {
    
    cerr << enspan_help;
    cerr << "Error: input- and output file is required.\n";
    exit(1);
  }
  if (nib_dir.empty())
  {
    
    cerr << enspan_help;
    cerr << "Error: nib file's root dir is required.\n";
    exit(1);
  }
  
  map<string, vector<discordant_pair>> enspan_map;
  string tmp;
  stringstream line;
  vector<double> insert;
  cout << "start to stats the insert size...\n";
  get_mean_insert_size(inp_file, insert);
  double insert_size = insert[0];
  double insert_sd = insert[1];
  cout << "the insert size mean: " << insert_size << ", the insert size sd:" << insert_sd << " .\n";
  double scan_dist, mask_dist, span_dist, cluster_dist;
  cluster_dist = span_dist = mask_dist = scan_dist = times * sqrt(times) * (insert_size + 3 * insert_sd);
  cout << "cluster_dist = span_dist = mask_dist = scan_dist = " << cluster_dist << " .\n";
  //cluster index
  vector<vector<long> > cluster_index;
  vector<cluster_info> cluster;
  vector<cluster_info> tmp_cluster_vec;
  int root_cluster_num;
  int removed_isolated_pair_count = 0;
  int after_cluster_count = 0;
  int scan_pairs_count = 0;
  scan_start = clock();
  scan_discordant_pairs(inp_file, build, qual, scan_dist, enspan_map, nib_dir);
//  cout << "discordant pairs found: " << enspan_vec.size() << endl;
  scan_end = clock();
  cout << "the scan discordant pairs step  costs time: " <<
       (scan_end - scan_start) / double(CLOCKS_PER_SEC) << " seconds \n";
  for (auto &chr_it : enspan_map)
  {
    cerr << "Now start to process the region: " << chr_it.first << endl;
    add_enspan_point_id(chr_it.second);
    remove_isolated_pairs(chr_it.second, mask_dist);
    
    if (chr_it.second.size() >= 2)
    {
      cluster_start = clock();
      removed_isolated_pair_count += (int) chr_it.second.size();
      if (fast_cluster)
      {
        root_cluster_num = find_cluster_pairs_enspan_fast(chr_it.second, cluster_dist, min_reads_per_cluster);
      }
      else
      {
        root_cluster_num = find_cluster_pairs_enspan_ahc(chr_it.second, cluster_dist, distance_type,
                                                         min_reads_per_cluster);
      }
      cout << "discordant pairs found after clustering: " << chr_it.second.size() << endl;
      cluster_end = clock();
      std::cout << "the current number of root cluster is:" << root_cluster_num << std::endl;
      cout << "the clustering  step  costs time: " <<
           (cluster_end - cluster_start) / double(CLOCKS_PER_SEC) << " seconds " << std::endl;
      //sort according to the cluster id
      sort(chr_it.second.begin(), chr_it.second.end(), cmp_enspan_id);
      //search for spanning reads
      bp_find_start = clock();
      vector<bam1_t *> split_reads;
      split_reads.clear();
      
      /// use SA tag to find breakpoint
      
      findClusterBreakPointInfoSaTag(inp_file, chr_it.second, span_dist, tmp_cluster_vec, split_reads, nib_dir);
      
      if (!tmp_cluster_vec.empty())
      {
        for (const auto &i : tmp_cluster_vec)
        {
          cluster.push_back(i);
        }
        tmp_cluster_vec.clear();
      }
      
    }
    cout << "discordant pairs found after removing isolated pairs: " << chr_it.second.size() << endl;
    
    
  }
//write enspan-file out
  write_enspan_out(out_file, cluster, filter);
  write_enspan_params(inp_file, out_file, build, cluster_dist, qual);
  end = clock();
  cout << "the fusion process of file " << inp_file << "  costs time: " <<
       (end - start) / double(CLOCKS_PER_SEC) << " seconds" << endl;
//  stat out
  ofstream out;
  // Append parameters to parameter file
  out.open((out_file + "_performance.txt").c_str());
  out
    << "scan_dist\tdiscordant pairs\tremove isolated\tafter_cluster\troot cluster\tscanning time\tcluster time\tfind breakpoint time\ttotal time"
    << std::endl;
  out << scan_dist << "\t";
  out << scan_pairs_count << "\t";
  out << removed_isolated_pair_count << "\t";
  out << after_cluster_count << "\t";
  out << root_cluster_num << "\t";
  out << (scan_end - scan_start) / double(CLOCKS_PER_SEC) << "\t";
  out << (cluster_end - cluster_start) / double(CLOCKS_PER_SEC) << "\t";
  out << (bp_find_end - bp_find_start) / double(CLOCKS_PER_SEC) << "\t";
  out << (end - start) / double(CLOCKS_PER_SEC);
  out << std::endl;
  out.close();
}

/**
 * @brief main function to make the cluster info
 * @param bam_file  input bam file
 * @param enspan  input spanning reads pairs vector
 * @param w half distance to extend to pull down the encompassing reads across the region
 * @param cluster_vec the produced cluster vector
 */
void
findClusterBreakPointInfoSaTag(string bam_file, vector<discordant_pair> &enspan, double w,
                               vector<cluster_info> &cluster_vec, vector<bam1_t *> &split_reads, string nib_dir)
{
  stringstream line;
  int k, current_index;
  cluster_info tmp_cluster;
  
  set<string> tmp_set;
  
  map<long, cluster_info> cluster_map;
  map<long, cluster_info>::iterator cluster_map_it;
  map<long, vector<int>> cluster_enspan_index_map;
  map<long, vector<int>>::iterator cluster_enspan_index_map_it;
  map<long, set<string>> cluster_enspan_type_map;
  map<long, set<string>>::iterator cluster_enspan_type_map_it;
  vector<int> enspan_index_vec;
  vector<int> tmp_vec(1, -1);
  int64_t mean_pos_dist = 0;
  
  
  //// collect index of enspan pairs for each cluster id
  if (!enspan.empty())
  {
    for (int k = 0; k < enspan.size(); k++)
    {
      cluster_enspan_type_map_it = cluster_enspan_type_map.find(enspan[k].cluster);
      if (cluster_enspan_type_map_it == cluster_enspan_type_map.end())
      {
        tmp_set.clear();
        if (enspan[k].p1_chr != enspan[k].p2_chr)
        {
          tmp_set.insert("diff_chr");
        }
        else
        {
          if (enspan[k].p1_strand == '-' and enspan[k].p2_strand == '+')
          {
            tmp_set.insert("same_chr_with_absolute_reverse");
            
          }
          
          if (enspan[k].p1_strand == enspan[k].p2_strand)
          {
            tmp_set.insert("same_chr_with_same_orientation");
            
          }
          
          if (enspan[k].p1_strand == '+' and enspan[k].p2_strand == '-')
          {
            tmp_set.insert("same_chr_with_default_orientation");
            
          }
        }
        
        cluster_enspan_type_map[enspan[k].cluster] = tmp_set;
      }
      else
      {
        if (enspan[k].p1_chr != enspan[k].p2_chr)
        {
          cluster_enspan_type_map[enspan[k].cluster].insert("diff_chr");
        }
        else
        {
          if (enspan[k].p1_strand == '-' and enspan[k].p2_strand == '+')
          {
            cluster_enspan_type_map[enspan[k].cluster].insert("same_chr_with_absolute_reverse");
            
          }
          
          if (enspan[k].p1_strand == enspan[k].p2_strand)
          {
            cluster_enspan_type_map[enspan[k].cluster].insert("same_chr_with_same_orientation");
            
          }
          
          if (enspan[k].p1_strand == '+' and enspan[k].p2_strand == '-')
          {
            cluster_enspan_type_map[enspan[k].cluster].insert("same_chr_with_default_orientation");
            
          }
        }
      }
      cluster_enspan_index_map_it = cluster_enspan_index_map.find(enspan[k].cluster);
      if (cluster_enspan_index_map_it == cluster_enspan_index_map.end())
      {
        tmp_vec[0] = k;
        cluster_enspan_index_map[enspan[k].cluster] = tmp_vec;
      }
      else
      {
        cluster_enspan_index_map[enspan[k].cluster].push_back(k);
      }
    }
    
    //// collect each cluster basic info: p1/2_chr, p1/2_mean_pos etc.
    for (cluster_enspan_index_map_it = cluster_enspan_index_map.begin();
      cluster_enspan_index_map_it != cluster_enspan_index_map.end(); ++cluster_enspan_index_map_it)
    {
      //// cluster initiation
      tmp_cluster.id = cluster_enspan_index_map_it->first;
      current_index = cluster_enspan_index_map_it->second[0];
      tmp_cluster.p1_chr = enspan[current_index].p1_chr;
      tmp_cluster.p2_chr = enspan[current_index].p2_chr;
      tmp_cluster.p1_mean_pos = enspan[current_index].p1_pos;
      tmp_cluster.p2_mean_pos = enspan[current_index].p2_pos;
      tmp_cluster.p1_min_pos = enspan[current_index].p1_pos;
      tmp_cluster.p2_min_pos = enspan[current_index].p2_pos;
      tmp_cluster.p1_max_pos = enspan[current_index].p1_pos;
      tmp_cluster.p2_max_pos = enspan[current_index].p2_pos;
      tmp_cluster.n_split_read = 0;
      tmp_cluster.n_discordant_pair = 1;
      tmp_cluster.discordant_reads = enspan[current_index].qname + ",";
      tmp_cluster.split_reads = "";
      tmp_cluster.p1_exact_pos = -1;
      tmp_cluster.p2_exact_pos = -1;
      tmp_cluster.inv = false;
      
      if (cluster_enspan_index_map_it->second.size() >= 2)
      {////only retain cluster that contain more than 2 spanning pairs.
        for (k = 1; k < cluster_enspan_index_map_it->second.size(); ++k)
        {
          current_index = cluster_enspan_index_map_it->second[k];
          tmp_cluster.p1_mean_pos += enspan[current_index].p1_pos;
          if (tmp_cluster.p1_min_pos > enspan[current_index].p1_pos)
            tmp_cluster.p1_min_pos = enspan[current_index].p1_pos;
          if (tmp_cluster.p1_max_pos < enspan[current_index].p1_pos)
            tmp_cluster.p1_max_pos = enspan[current_index].p1_pos;
          
          tmp_cluster.p2_mean_pos += enspan[current_index].p2_pos;
          if (tmp_cluster.p2_min_pos > enspan[current_index].p2_pos)
            tmp_cluster.p2_min_pos = enspan[current_index].p2_pos;
          if (tmp_cluster.p2_max_pos < enspan[current_index].p2_pos)
            tmp_cluster.p2_max_pos = enspan[current_index].p2_pos;
          
          tmp_cluster.discordant_reads += enspan[current_index].qname + ",";
          tmp_cluster.n_discordant_pair++;
        }
      }
      
      tmp_cluster.p1_mean_pos = (uint32_t) ((double) tmp_cluster.p1_mean_pos / (double) tmp_cluster.n_discordant_pair);
      tmp_cluster.p2_mean_pos = (uint32_t) ((double) tmp_cluster.p2_mean_pos / (double) tmp_cluster.n_discordant_pair);
      
      mean_pos_dist = (int64_t) (tmp_cluster.p1_mean_pos - tmp_cluster.p2_mean_pos);
      tmp_cluster.drp_type_set = cluster_enspan_type_map[tmp_cluster.id];
      
      if (!(tmp_cluster.p1_chr == tmp_cluster.p2_chr && mean_pos_dist <= 2 * w && mean_pos_dist >= -2 * w))
      {
        cluster_map[tmp_cluster.id] = tmp_cluster; ////only retain cluster whose p1_pos and p2_pos is not so close.
      }
    }
    cout << "there is " << cluster_map.size() << " cluster after produce cluster data\n";
    //// determine breakpoint pos for each cluster, the result is recorded in the object cluster map
    findEncompassingReadsAndBreakPointInfo(cluster_map, bam_file, w, split_reads);
    cluster_vec.clear();
    for (cluster_map_it = cluster_map.begin(); cluster_map_it != cluster_map.end(); ++cluster_map_it)
    {
      cluster_vec.push_back(cluster_map_it->second); //// convert cluster map to cluster vector
    }
    
    map<long, vector<discordant_pair>> cluster_enspan_vec_map;
    map<long, vector<discordant_pair>>::iterator cluster_enspan_vec_map_it;
    vector<discordant_pair> tmp_pair;
    int enspan_index;
    for (cluster_enspan_index_map_it = cluster_enspan_index_map.begin();
      cluster_enspan_index_map_it != cluster_enspan_index_map.end(); ++cluster_enspan_index_map_it)
    {
      tmp_pair.clear();
      for (int i = 0; i < cluster_enspan_index_map_it->second.size(); ++i)
      {
        enspan_index = cluster_enspan_index_map_it->second[i];
        tmp_pair.push_back(enspan[enspan_index]);
      }
      cluster_enspan_vec_map[cluster_enspan_index_map_it->first] = tmp_pair;
    }
    
    annotate_cluster_for_sa_tag(cluster_vec, nib_dir);
    
  }
}

/**
 * @brief the module that determine the breakpoint pos for both part of each fusion cluster
 * @param cluster_map: the cluster map which will record the info of breakpoint poses
 * @param bam_file: the input bam file
 * @param w : half of the region span that will be pull down the encompassing reads(spilt reads that across breakpoint)
 */
void findEncompassingReadsAndBreakPointInfo(map<long, cluster_info> &cluster_map, const string bam_file, const int w,
                                            vector<bam1_t *> &split_reads)
{
  uint32_t p1_region_start, p1_region_end, p2_region_start, p2_region_end;
  string p1_chr, p2_chr;
  map<long, cluster_info>::iterator cluster_map_it;
  map<long, cluster_info>::iterator cluster_map_new_it;
  map<long, cluster_info> tmp_cluster_map;
  bam_index_t *idx_bam;
  map<string, vector<split_align_pair>> p1_encompass_map, p2_encompass_map;
  breakpoint_pair bp_pair;
  
  /**
   * open the bam file
   */
  samfile_t *fp = samopen(bam_file.c_str(), "rb", 0);
  if (fp == NULL)
  {
    std::cerr << "Error: can not open input bam-file:\t" << bam_file << std::endl;
    exit(1);
  }
  idx_bam = bam_index_load(bam_file.c_str());
  if (idx_bam == NULL)
  {
    std::cerr << "Error: please index bam-file first:\t" << bam_file << std::endl;
    exit(1);
  }
  int valid_count = 0;
  /**
   * iter the cluster map, detect the real info
   */
  for (cluster_map_it = cluster_map.begin(); cluster_map_it != cluster_map.end(); ++cluster_map_it)
  {
    p1_encompass_map.clear();
    p2_encompass_map.clear();
    bp_pair.encompass_num = 0;
    bp_pair.p1_chr = "";
    bp_pair.p1_bp = -1;
    bp_pair.p2_chr = "";
    bp_pair.p2_bp = -1;
    p1_region_start = (uint32_t) (cluster_map_it->second.p1_mean_pos - w);
    p1_region_end = (uint32_t) (cluster_map_it->second.p1_mean_pos + w);
    p2_region_start = (uint32_t) (cluster_map_it->second.p2_mean_pos - w);
    p2_region_end = (uint32_t) (cluster_map_it->second.p2_mean_pos + w);
    p1_chr = cluster_map_it->second.p1_chr;
    p2_chr = cluster_map_it->second.p2_chr;
    find_sa_reads(fp, p1_chr, p1_region_start, p1_region_end, p1_encompass_map, idx_bam);
    /* if split reads found in p1 part, then find split reads in p2 part*/
    if (p1_encompass_map.size() > 0)
      find_sa_reads(fp, p2_chr, p2_region_start, p2_region_end, p2_encompass_map, idx_bam);
    
    ////only retain the cluster that both side has split reads
    if (p1_encompass_map.size() > 0 && p2_encompass_map.size() > 0)
    {
      find_bp_pair(p1_encompass_map, p2_encompass_map, bp_pair, p1_chr, p2_chr, split_reads,
                   2);//dig the breakpoint info from encompass map in both side
      if (bp_pair.encompass_num >= 2)
      {
        cluster_map_it->second.p1_exact_pos = bp_pair.p1_bp;
        cluster_map_it->second.p2_exact_pos = bp_pair.p2_bp;
        cluster_map_it->second.n_split_read = bp_pair.encompass_num;
        cluster_map_it->second.p1_coverage = cal_mean_depth_oc(cluster_map_it->second.p1_chr,
                                                               cluster_map_it->second.p1_min_pos >
                                                               cluster_map_it->second.p1_exact_pos
                                                               ? cluster_map_it->second.p1_exact_pos
                                                               : cluster_map_it->second.p1_min_pos,
                                                               cluster_map_it->second.p1_max_pos >
                                                               cluster_map_it->second.p1_exact_pos
                                                               ? cluster_map_it->second.p1_max_pos
                                                               : cluster_map_it->second.p1_exact_pos, fp,
                                                               idx_bam);
        cluster_map_it->second.p2_coverage = cal_mean_depth_oc(cluster_map_it->second.p2_chr,
                                                               cluster_map_it->second.p2_min_pos >
                                                               cluster_map_it->second.p2_exact_pos
                                                               ? cluster_map_it->second.p2_exact_pos
                                                               : cluster_map_it->second.p2_min_pos,
                                                               cluster_map_it->second.p2_max_pos >
                                                               cluster_map_it->second.p2_exact_pos
                                                               ? cluster_map_it->second.p2_max_pos
                                                               : cluster_map_it->second.p2_exact_pos, fp,
                                                               idx_bam);
        cluster_map_it->second.p1_bp_depth = cal_single_base_depth(cluster_map_it->second.p1_chr,
                                                                   cluster_map_it->second.p1_exact_pos, fp, idx_bam);
        cluster_map_it->second.p2_bp_depth = cal_single_base_depth(cluster_map_it->second.p2_chr,
                                                                   cluster_map_it->second.p2_exact_pos, fp, idx_bam);
        cluster_map_it->second.p1_alle_freq =
          (float) cluster_map_it->second.n_split_read / (float) cluster_map_it->second.p1_bp_depth;
        cluster_map_it->second.p2_alle_freq =
          (float) cluster_map_it->second.n_split_read / (float) cluster_map_it->second.p2_bp_depth;
        cluster_map_it->second.fusion_type = determine_fusion_type_from_drp(cluster_map_it->second);
        
        tmp_cluster_map[cluster_map_it->first] = cluster_map_it->second;
        valid_count++;
      }
    }
  }
  cluster_map = tmp_cluster_map;
  tmp_cluster_map.clear();
  cout << "valid cluster count: " << valid_count << endl;
  samclose(fp);
}

void annotate_cluster_for_sa_tag(vector<cluster_info> &clusters, string nib_dir)
{
  string ref_gene = (std::string) INSTALLDIR + "/ref_files/refGene.txt";  /// refGene.txt
  map<string, int> hotspot_fusions;
  map<string, int>::iterator hotspot_fusions_it;
  vector<map<string, string>> sino_anno_fusion_info;
  vector<map<string, string>> cosmic_anno_fusion_info;
  
  std::vector<RefSeqTranscript> txpts;
  readRefSeqTranscript(ref_gene, txpts);/// get all transcripts from ref gene file
  std::vector<RefSeqTranscript> longest_txpts;
  vector<string> gene_vec;
  string fusion_5_3_exon_pair, fusion_5exon_pair, fusion_3exon_pair, fusion_no_exon_pair;
  
  add_cds_parts(txpts);/// add cds info to all txpts
  long p1_pos, p2_pos;
  string chr1, chr2;
  
  vector<string> exon_infos;
  vector<cluster_info> temp_clusters;
  
  for (auto cluster:clusters)
  {
    chr1 = cluster.p1_chr;
    chr2 = cluster.p2_chr;
    
    if (cluster.p1_exact_pos == -1)
    {
      p1_pos = cluster.p1_mean_pos;
    }
    else
    {
      p1_pos = cluster.p1_exact_pos;
    }
    
    if (cluster.p2_exact_pos == -1)
    {
      p2_pos = cluster.p2_mean_pos;
    }
    else
    {
      p2_pos = cluster.p2_exact_pos;
    }
    add_exon_anno(txpts, chr1, chr2, p1_pos, p2_pos,
                  exon_infos, cluster.p1_part, cluster.p2_part);
    /// load the return value
    cluster.p1_behalf_gene = exon_infos[0];
    cluster.p2_behalf_gene = exon_infos[1];
    cluster.p1_exon_info = exon_infos[2];
    cluster.p2_exon_info = exon_infos[3];
    cluster.p1_strand = exon_infos[4];
    cluster.p2_strand = exon_infos[5];
    cluster.p1_genes = exon_infos[6];
    cluster.p2_genes = exon_infos[7];
    cluster.p1_gene_part = exon_infos[8];
    cluster.p2_gene_part = exon_infos[9];
    cluster.p1_bp_exon = exon_infos[10];
    cluster.p2_bp_exon = exon_infos[11];
    cluster.up_gene = exon_infos[12];
    cluster.down_gene = exon_infos[13];
    cluster.fusion_pair = exon_infos[14];
    
    string p1_left_seq = get_left_neighbor_sequence_nib(cluster.p1_chr, cluster.p1_exact_pos, 20, nib_dir);
    string p1_right_seq = get_right_neighbor_sequence_nib(cluster.p1_chr, cluster.p1_exact_pos - 1, 21, nib_dir);
    string p2_left_seq = get_left_neighbor_sequence_nib(cluster.p2_chr, cluster.p2_exact_pos, 20, nib_dir);
    string p2_right_seq = get_right_neighbor_sequence_nib(cluster.p2_chr, cluster.p2_exact_pos - 1, 21, nib_dir);
    cluster.p1_rpt = p1_left_seq + p1_right_seq;
    cluster.p2_rpt = p2_left_seq + p2_right_seq;
    cluster.is_rpt = (find_longest_repeat_substring(cluster.p1_rpt) > 10 ||
                      find_longest_repeat_substring(cluster.p2_rpt) > 10);
    temp_clusters.push_back(cluster);
    exon_infos.clear();
  }
  clusters = temp_clusters;
  temp_clusters.clear();
}

/**
 * @brief dig breakpoint poses of the current cluster from split align pair map
 * @param p1_encompass_map:the encompass map of p1 part
 * @param p2_encompass_map:the encompass map of p2 part
 * @param bp_pair:the produced breakpoint info struct
 * @param p1_chr:the chromosome of p1 part
 * @param p2_chr:the chromosome of p2 part
 */
void find_bp_pair(map<string, vector<split_align_pair>> &p1_encompass_map,
                  map<string, vector<split_align_pair>> &p2_encompass_map, breakpoint_pair &bp_pair,
                  const string &p1_chr, const string &p2_chr, vector<bam1_t *> &split_reads, int bp_pos_error)
{
  bam1_t *b_copy;
  map<string, vector<split_align_pair>>::iterator p1_encompass_map_it, p2_encompass_map_it;
  bp_pair.p1_chr = p1_chr;
  bp_pair.p2_chr = p2_chr;
  split_align_pair p1_tmp_align_pair, p2_tmp_align_pair;
  map<uint32_t, int> p1_bp_pos_count_map, p2_bp_pos_count_map;
  map<uint32_t, int>::iterator p1_bp_pos_count_it, p2_bp_pos_count_it;
  vector<string> sa_substrs;
  bool condition, new_condition;
 // int max_count, max_count_new;
  vector<pair<int32_t, int32_t >> bp_pos_pair_vec;
  vector<bp_pos_pair> bp_pos_pair_vec_new;
  vector<bp_pos_pair> bp_pos_pair_vec_update;
  bp_pos_pair pair1, pair2, pair3, pair4, pair_update;
  bp_pos_pair_vec.clear();
  p1_bp_pos_count_map.clear();
  p2_bp_pos_count_map.clear();
  map<string, int> bp_pos_pair_count_map, bp_part_pair_count_map;
  map<string, int> bp_pos_pair_count_map_update, bp_part_pair_count_map_update, bp_pos_pair_count_map_update_temp;
  map<string, int>::iterator string_int_map_it;
  for (p1_encompass_map_it = p1_encompass_map.begin();
    p1_encompass_map_it != p1_encompass_map.end(); ++p1_encompass_map_it)
  {
    p2_encompass_map_it = p2_encompass_map.find(p1_encompass_map_it->first);
    if (p2_encompass_map_it == p2_encompass_map.end())
      continue;
    else
      ////only consider the reads existed in both side
    {
      for (int i = 0; i < p1_encompass_map_it->second.size(); ++i)
      {
        for (int j = 0; j < p2_encompass_map_it->second.size(); ++j)
        {
          p1_tmp_align_pair = p1_encompass_map_it->second[i];
          p2_tmp_align_pair = p2_encompass_map_it->second[j];
          condition = (p1_tmp_align_pair.secondary != p2_tmp_align_pair.secondary) &&
                      (p1_tmp_align_pair.primary_chr == p2_tmp_align_pair.primary_chr) &&
                      (p1_tmp_align_pair.secondary_chr == p2_tmp_align_pair.secondary_chr) &&
                      (p1_tmp_align_pair.primary_start == p2_tmp_align_pair.primary_start) &&
                      (p1_tmp_align_pair.secondary_start == p2_tmp_align_pair.secondary_start) &&
                      (p1_tmp_align_pair.primary_end == p2_tmp_align_pair.primary_end) &&
                      (p1_tmp_align_pair.secondary_end == p2_tmp_align_pair.secondary_end) &&
                      (p1_tmp_align_pair.primary_cigar_str == p2_tmp_align_pair.primary_cigar_str) &&
                      (p1_tmp_align_pair.secondary_cigar_str == p2_tmp_align_pair.secondary_cigar_str);
          
          new_condition = (p1_tmp_align_pair.secondary != p2_tmp_align_pair.secondary) &&
                          (p1_tmp_align_pair.primary_chr == p2_tmp_align_pair.primary_chr) &&
                          (p1_tmp_align_pair.secondary_chr == p2_tmp_align_pair.secondary_chr) &&
                          (p1_tmp_align_pair.primary_start == p2_tmp_align_pair.primary_start) &&
                          (p1_tmp_align_pair.secondary_start == p2_tmp_align_pair.secondary_start) &&
                          (p1_tmp_align_pair.primary_end == p2_tmp_align_pair.primary_end) &&
                          (p1_tmp_align_pair.secondary_end == p2_tmp_align_pair.secondary_end) &&
                          (p1_tmp_align_pair.primary_cigar_str == p2_tmp_align_pair.primary_cigar_str) &&
                          (p1_tmp_align_pair.secondary_cigar_str == p2_tmp_align_pair.secondary_cigar_str) &&
                          (p1_tmp_align_pair.primary_bp == p2_tmp_align_pair.primary_bp) &&
                          (p1_tmp_align_pair.secondary_bp == p2_tmp_align_pair.secondary_bp);
          
          if (new_condition)
          {
            b_copy = bam_init1();
            split_reads.push_back(bam_copy1(b_copy, p1_tmp_align_pair.current_align));
            
            b_copy = bam_init1();
            split_reads.push_back(bam_copy1(b_copy, p2_tmp_align_pair.current_align));
            
            if (p1_tmp_align_pair.primary_chr == p1_chr)
              ////primary part of the reads is in p1 part
            {
              
              pair1.p1_pos = (int32_t) p1_tmp_align_pair.primary_start;
              pair1.p2_pos = (int32_t) p1_tmp_align_pair.secondary_start;
              pair1.p1_match_part = "right", pair1.p2_match_part = "right";
              
              pair2.p1_pos = (int32_t) p1_tmp_align_pair.primary_start;
              pair2.p2_pos = (int32_t) p1_tmp_align_pair.secondary_end;
              pair2.p1_match_part = "right", pair2.p2_match_part = "left";
              
              pair3.p1_pos = (int32_t) p1_tmp_align_pair.primary_end;
              pair3.p2_pos = (int32_t) p1_tmp_align_pair.secondary_start;
              pair3.p1_match_part = "left", pair3.p2_match_part = "right";
              
              pair4.p1_pos = (int32_t) p1_tmp_align_pair.primary_end;
              pair4.p2_pos = (int32_t) p1_tmp_align_pair.secondary_end;
              pair4.p1_match_part = "left", pair4.p2_match_part = "left";
              bp_pos_pair_vec_new.push_back(pair1);
              bp_pos_pair_vec_new.push_back(pair2);
              bp_pos_pair_vec_new.push_back(pair3);
              bp_pos_pair_vec_new.push_back(pair4);
              
              pair_update.p1_pos = (int32_t) p1_tmp_align_pair.primary_bp;
              pair_update.p2_pos = (int32_t) p1_tmp_align_pair.secondary_bp;
              
              if (p1_tmp_align_pair.primary_bp == p1_tmp_align_pair.primary_start)
                pair_update.p1_match_part = "right";
              else
                pair_update.p1_match_part = "left";
              
              if (p1_tmp_align_pair.secondary_bp == p1_tmp_align_pair.secondary_start)
                pair_update.p2_match_part = "right";
              else
                pair_update.p2_match_part = "left";
              bp_pos_pair_vec_update.push_back(pair_update);
            }
            else
            {
              ////primary part of the reads is in p2 part
//              bp_pos_pair_vec.push_back(
//                make_pair((int32_t) p1_tmp_align_pair.secondary_start, (int32_t) p1_tmp_align_pair.primary_start));
//              bp_pos_pair_vec.push_back(
//                make_pair((int32_t) p1_tmp_align_pair.secondary_start, (int32_t) p1_tmp_align_pair.primary_end));
//              bp_pos_pair_vec.push_back(
//                make_pair((int32_t) p1_tmp_align_pair.secondary_end, (int32_t) p1_tmp_align_pair.primary_start));
//              bp_pos_pair_vec.push_back(
//                make_pair((int32_t) p1_tmp_align_pair.secondary_end, (int32_t) p1_tmp_align_pair.primary_end));
              
              pair1.p1_pos = (int32_t) p1_tmp_align_pair.secondary_start;
              pair1.p2_pos = (int32_t) p1_tmp_align_pair.primary_start;
              pair1.p1_match_part = "right", pair1.p2_match_part = "right";
              
              pair2.p1_pos = (int32_t) p1_tmp_align_pair.secondary_start;
              pair2.p2_pos = (int32_t) p1_tmp_align_pair.primary_end;
              pair2.p1_match_part = "right", pair2.p2_match_part = "left";
              
              pair3.p1_pos = (int32_t) p1_tmp_align_pair.secondary_end;
              pair3.p2_pos = (int32_t) p1_tmp_align_pair.primary_start;
              pair3.p1_match_part = "left", pair3.p2_match_part = "right";
              
              pair4.p1_pos = (int32_t) p1_tmp_align_pair.secondary_end;
              pair4.p2_pos = (int32_t) p1_tmp_align_pair.primary_end;
              pair4.p1_match_part = "left", pair4.p2_match_part = "left";
              bp_pos_pair_vec_new.push_back(pair1);
              bp_pos_pair_vec_new.push_back(pair2);
              bp_pos_pair_vec_new.push_back(pair3);
              bp_pos_pair_vec_new.push_back(pair4);
              
              pair_update.p1_pos = (int32_t) p1_tmp_align_pair.secondary_bp;
              pair_update.p2_pos = (int32_t) p1_tmp_align_pair.primary_bp;

//              if (p1_tmp_align_pair.secondary_bp==p1_tmp_align_pair.primary_start)
//                pair_update.p1_match_part = "right";
//              else
//                pair_update.p1_match_part = "left";
//
//              if (p1_tmp_align_pair.primary_bp==p1_tmp_align_pair.secondary_start)
//                pair_update.p2_match_part = "right";
//              else
//                pair_update.p2_match_part = "left";
//              pair_update.p1_match_part = "right";
//              pair_update.p1_match_part = "";
              
              bp_pos_pair_vec_update.push_back(pair_update);
            }
          }
        }
      }
    }
  }
  ///////////////////////////////////////////////////////////old version///////////////////////////////////////
/**
  max_count_new = 0;
  string tmp_str_new = "";
  
  vector<string> tmp_substrs_new;
  for (int k = 0; k < bp_pos_pair_vec_new.size(); ++k)
  {
    tmp_str_new = to_string(bp_pos_pair_vec_new[k].p1_pos) + "," + to_string(bp_pos_pair_vec_new[k].p2_pos);
    string_int_map_it = bp_pos_pair_count_map.find(tmp_str_new);
    if (string_int_map_it == bp_pos_pair_count_map.end())
      bp_pos_pair_count_map[tmp_str_new] = 1;
    else
      string_int_map_it->second++;
  }
  for (string_int_map_it = bp_pos_pair_count_map.begin();
    string_int_map_it !=
    bp_pos_pair_count_map.end(); ++string_int_map_it)/// stat the repeat count of each breakpoint pos pair
  {
    if (max_count_new < string_int_map_it->second)
    {
      max_count_new = string_int_map_it->second;
      tmp_substrs_new = split_string(string_int_map_it->first, ",");
      bp_pair.p1_bp = (uint32_t) stoull(tmp_substrs_new[0]);
      bp_pair.p2_bp = (uint32_t) stoull(tmp_substrs_new[1]);
    }
  }
  bp_pair.encompass_num = max_count_new;
  
  bp_part_pair_count_map.clear();
  for (int m = 0; m < bp_pos_pair_vec_new.size(); ++m)
  {
    if ((bp_pair.p1_bp == bp_pos_pair_vec_new[m].p1_pos) && (bp_pair.p2_bp == bp_pos_pair_vec_new[m].p2_pos))
    {
      tmp_str_new = bp_pos_pair_vec_new[m].p1_match_part + "," + bp_pos_pair_vec_new[m].p2_match_part;
      string_int_map_it = bp_part_pair_count_map.find(tmp_str_new);
      if (string_int_map_it == bp_part_pair_count_map.end())
        bp_part_pair_count_map[tmp_str_new] = 1;
      else
        string_int_map_it->second++;
    }
  }
  max_count_new = 0;
  for (string_int_map_it = bp_part_pair_count_map.begin();
    string_int_map_it != bp_part_pair_count_map.end(); ++string_int_map_it)
  {
    if (max_count_new < string_int_map_it->second)
    {
      max_count_new = string_int_map_it->second;
      tmp_substrs_new = split_string(string_int_map_it->first, ",");
      bp_pair.p1_part = tmp_substrs_new[0];
      bp_pair.p2_part = tmp_substrs_new[1];
    }
  }
  */

///////////////////////////////////////////////////////////////////////update  version///////////////////////////////
  
  int max_count_update = 0;
  string tmp_str_update = "";
  vector<string> tmp_substrs_update;
  uint32_t tmp_p1_pos, tmp_p2_pos;
  bp_pos_pair tmp_bps;
  
  
  for (int k = 0; k < bp_pos_pair_vec_update.size(); ++k)
  {
    tmp_str_update = to_string(bp_pos_pair_vec_update[k].p1_pos) + "," + to_string(bp_pos_pair_vec_update[k].p2_pos);
    bp_pos_pair_count_map_update[tmp_str_update] = 0;
  }
  
  for (string_int_map_it = bp_pos_pair_count_map_update.begin();
    string_int_map_it !=
    bp_pos_pair_count_map_update.end(); ++string_int_map_it)
  {
    tmp_substrs_update = split_string(string_int_map_it->first, ",");
    tmp_p1_pos = (uint32_t) stoull(tmp_substrs_update[0]);
    tmp_p2_pos = (uint32_t) stoull(tmp_substrs_update[1]);
    for (int k = 0; k < bp_pos_pair_vec_update.size(); ++k)
    {
      tmp_bps = bp_pos_pair_vec_update[k];
      if ((tmp_bps.p1_pos <= tmp_p1_pos + bp_pos_error && tmp_bps.p1_pos >= tmp_p1_pos - bp_pos_error) &&
          (tmp_bps.p2_pos <= tmp_p2_pos + bp_pos_error && tmp_bps.p2_pos >= tmp_p2_pos - bp_pos_error))
      {
        string_int_map_it->second++;
      }
    }
  }


//  for (int k = 0; k < bp_pos_pair_vec_update.size(); ++k)
//  {
//    tmp_str_update = to_string(bp_pos_pair_vec_update[k].p1_pos) + "," + to_string(bp_pos_pair_vec_update[k].p2_pos);
//
//    string_int_map_it = bp_pos_pair_count_map_update.find(tmp_str_update);
//    if (string_int_map_it == bp_pos_pair_count_map_update.end())
//      bp_pos_pair_count_map_update[tmp_str_update] = 1;
//    else
//      string_int_map_it->second++;
//  }
//
  
  for (string_int_map_it = bp_pos_pair_count_map_update.begin();
    string_int_map_it !=
    bp_pos_pair_count_map_update.end(); ++string_int_map_it)/// stat the repeat count of each breakpoint pos pair
  {
    if (max_count_update < string_int_map_it->second)
    {
      max_count_update = string_int_map_it->second;
      tmp_substrs_update = split_string(string_int_map_it->first, ",");
      bp_pair.p1_bp = (uint32_t) stoull(tmp_substrs_update[0]);
      bp_pair.p2_bp = (uint32_t) stoull(tmp_substrs_update[1]);
//      bp_pair.p1_part = tmp_substrs_update[2];
//      bp_pair.p2_part = tmp_substrs_update[3];
    }
  }
  bp_pair.encompass_num = max_count_update;
  
}

/**
 * @brief pull down all the split reads alignment, and extract split align pair info from the select alignment
 * @param fp : bam file reader
 * @param region_chr : region chromosome
 * @param region_start : region start pos
 * @param region_end: region end pos
 * @param encompassing_map: the produced encompassing map
 * @param idx_bam : bam index
 */
void find_sa_reads(samfile_t *fp, const string region_chr, uint32_t region_start, uint32_t region_end,
                   map<string, vector<split_align_pair>> &encompassing_map, bam_index_t *idx_bam)
{
  int tid_t;
  int beg_t, end_t;
  bam_iter_t iter_bam;
  string sa_tag, oc_tag, read_name;
  BamAlignment align;
  uint32_t sa_start, sa_end;
  bool condition, advanced_contition;
  vector<string> sa_tag_substr;
  CigarRoller tmp_sa_cigar, tmp_cigar;
  bam_parse_region(fp->header, region_chr.c_str(), &tid_t, &beg_t, &end_t);
  iter_bam = bam_iter_query(idx_bam, tid_t, region_start, region_end);//// bam index
  bam1_t *b = bam_init1();
  int total_coverage = 0;
  int total_evidence_alignments_num = 0;
  map<string, vector<split_align_pair>>::iterator encompassing_map_it;
  encompassing_map.clear();
  vector<split_align_pair> tmp_split_align_pair_vector;
  split_align_pair tmp_pair;
  ///when determine the complementary cigar, we can tolerate mismatch numbers
  int mismatch_num = 10;
  while (bam_iter_read(fp->x.bam, iter_bam, b) >= 0)
  {
    total_coverage++;
    align.readBamRecord(b);
    sa_tag = align.saTag();
    //// only consider the reads that has sa tag
    condition = (sa_tag != "") && (b->core.qual >= 0) && (!(b->core.flag & BAM_FDUP)) && (b->core.flag & BAM_FPAIRED);
    
    if (condition)
    {
      oc_tag = align.originCigar();
      sa_tag_substr.clear();
      sa_tag_substr = split_string(sa_tag, ",");
      tmp_sa_cigar.Set(sa_tag_substr[3].c_str());
      if (!oc_tag.empty())
      {
        tmp_cigar.Set(oc_tag.c_str());
      }
      else
      {
        tmp_cigar.Set(align.getCigarString().c_str());
      }
      ///
      advanced_contition = tmp_cigar.is_complementary_cigar(sa_tag_substr[3], mismatch_num);
      
      if (advanced_contition)
      {
        total_evidence_alignments_num++;
        read_name = align.getReadName();
        encompassing_map_it = encompassing_map.find(read_name);
        tmp_pair.flag = b->core.flag;
        tmp_pair.read_name = read_name;
        tmp_pair.sa_tag = sa_tag;
        tmp_pair.secondary = (bool) (tmp_pair.flag & BAM_FSECONDARY); //// judge it is secondary alignment or not
        tmp_pair.current_align = bam_init1();
        tmp_pair.current_align = bam_copy1(tmp_pair.current_align, b);
        
        sa_start = (uint32_t) stoi(sa_tag_substr[1]);
        sa_end = (uint32_t) tmp_sa_cigar.getAlignmentEnd(sa_start);
        
        
        if (!tmp_pair.secondary)
        {
          tmp_pair.primary_chr = align.getChrName();
          tmp_pair.primary_start = (uint32_t) align.getAlignmentStart();
          
          if (oc_tag != "")
          {
            tmp_pair.primary_cigar_str = oc_tag;
            tmp_cigar.Set(oc_tag.c_str());
            tmp_pair.primary_end = (uint32_t) tmp_cigar.getAlignmentEnd(
              (uint32_t) align.getAlignmentStart());
          }
          else
          {
            tmp_pair.primary_cigar_str = align.getCigarString();
            tmp_pair.primary_end = (uint32_t) align.getAlignmentEnd();
          }
          if (tmp_cigar.getNumBeginClips() != 0)
            tmp_pair.primary_bp = (uint32_t) align.getAlignmentStart();
          else if (tmp_cigar.getNumEndClips() != 0)
            tmp_pair.primary_bp = (uint32_t) align.getAlignmentEnd();
          else
          {
            cerr << "error cigar: " << std::endl;
            exit(-1);
          }
          
          if (tmp_sa_cigar.getNumBeginClips() != 0)
            tmp_pair.secondary_bp = sa_start;
          else if (tmp_sa_cigar.getNumEndClips() != 0)
            tmp_pair.secondary_bp = sa_end;
          else
          {
            cerr << "error cigar: " << std::endl;
            exit(-1);
          }
          tmp_pair.secondary_chr = sa_tag_substr[0];
          tmp_pair.secondary_start = sa_start;
          tmp_pair.secondary_end = sa_end;
          tmp_pair.secondary_cigar_str = sa_tag_substr[3];
        }
        else
        {
          tmp_pair.primary_chr = sa_tag_substr[0];
          tmp_pair.primary_start = sa_start;
          tmp_pair.primary_end = sa_end;
          tmp_pair.primary_cigar_str = sa_tag_substr[3];
          if (tmp_sa_cigar.getNumBeginClips() != 0)
            tmp_pair.primary_bp = sa_start;
          else if (tmp_sa_cigar.getNumEndClips() != 0)
            tmp_pair.primary_bp = sa_end;
          else
          {
            cerr << "error cigar: " << std::endl;
            exit(-1);
          }
          
          if (oc_tag != "")
          {
            tmp_pair.secondary_cigar_str = oc_tag;
            tmp_cigar.Set(oc_tag.c_str());
            tmp_pair.secondary_end = (uint32_t) tmp_cigar.getAlignmentEnd(
              (uint32_t) align.getAlignmentStart());
          }
          else
          {
            tmp_pair.secondary_cigar_str = align.getCigarString();
            tmp_pair.secondary_end = (uint32_t) align.getAlignmentEnd();
          }
          
          if (tmp_cigar.getNumBeginClips() != 0)
            tmp_pair.secondary_bp = (uint32_t) align.getAlignmentStart();
          else if (tmp_cigar.getNumEndClips() != 0)
            tmp_pair.secondary_bp = (uint32_t) align.getAlignmentEnd();
          else
          {
            cerr << "error cigar: " << std::endl;
            exit(-1);
          }
          
          
          tmp_pair.secondary_chr = align.getChrName();
          tmp_pair.secondary_start = (uint32_t) align.getAlignmentStart();
        }
        
        if (encompassing_map_it == encompassing_map.end())
        {
          tmp_split_align_pair_vector.clear();
          tmp_split_align_pair_vector.push_back(tmp_pair); //// load the new split align pair
          encompassing_map[read_name] = tmp_split_align_pair_vector;
        }
        else
          encompassing_map_it->second.push_back(tmp_pair);
        
      }
      
    }
  }
  //// only consider the region that has minimum coverage of 5 and has minimum total split reads alignment of 2
  if (total_coverage < 5 || total_evidence_alignments_num < 2)
  {
    encompassing_map.clear();
  }
  return;
}

/**
 * fast cluster method that using sort enspan pairs vector
 * @param enspan : input discordant pair vector
 * @param w
 * @param min_reads : minimum pairs number of each cluster should have
 * @return
 */
int find_cluster_pairs_enspan_fast(vector<discordant_pair> &enspan, double w, int min_reads)
{
  int i, j, k, n;
  vector<int> cl_index;
  vector<discordant_pair> tmp_enspan;
  stringstream line;
  long pre_pos;
  map<string, int> key, key_cl;
  map<string, int>::iterator it, it2;
  //find clusters in first pair
  k = 1;
  n = (int) enspan.size();
  pre_pos = enspan[0].p1_chr_pos;
  cl_index.push_back(0);////cl_index: it contain the info :   enspan_pair indexes set which the current cluster has
  
  for (i = 1; i < n; i++)
  {
    ////check if reads of first pair are in window
    if (enspan[i].p1_chr_pos <= pre_pos + w && i != n - 1)
    {
      cl_index.push_back(i);
    }
    else
    {
      ////only report clusters that have more than 2 members
      if (cl_index.size() >= min_reads)
      {
        for (j = 0; j < cl_index.size(); j++)
        {
          line.str("");
          line.clear();
          line << k << ":"; ////cluster_id is k
          enspan[cl_index[j]].cluster_id = line.str();
          tmp_enspan.push_back(enspan[cl_index[j]]); ////only retain the pair which is included in a cluster
        }
        k++;
      }
      pre_pos = enspan[i].p1_chr_pos; //// change the  pre_pos to current pair's pos
      cl_index.clear();
      cl_index.push_back(i);
    }
  }
  enspan = tmp_enspan;
  tmp_enspan.clear();
  cl_index.clear();
  sort(enspan.begin(), enspan.end(), cmp_p2_enspan_pairs);
  //find clusters in second pair
  k = 1;
  n = (int) enspan.size();
  pre_pos = enspan[0].p2_chr_pos;
  cl_index.push_back(0);
  for (i = 1; i < n; i++)
  {
    //check if reads of first pair are in window
    if (enspan[i].p2_chr_pos <= pre_pos + w && i != n - 1)
    {
      cl_index.push_back(i);
    }
    else
    {
      //only report clusters that have more than 2 members
      if (cl_index.size() >= min_reads)
      {
        for (j = 0; j < cl_index.size(); j++)
        {
          line.str("");
          line.clear();
          line << enspan[cl_index[j]].cluster_id << k;
          enspan[cl_index[j]].cluster_id = line.str();
          tmp_enspan.push_back(enspan[cl_index[j]]);
        }
        k++;
      }
      pre_pos = enspan[i].p2_chr_pos;
      cl_index.clear();
      cl_index.push_back(i);
    }
  }
  enspan = tmp_enspan;
  tmp_enspan.clear();
  cl_index.clear();
  sort(enspan.begin(), enspan.end(), cmp_p1_enspan_pairs);
  //find matching cluster id
  for (i = 0; i < enspan.size(); i++)
  {
    it = key.find(enspan[i].cluster_id); //// key is record the pair count for each cluster_id
    if (it == key.end())
      key[enspan[i].cluster_id] = 1;
    else
      it->second++;
  }
  k = 0;
  key_cl.clear();
  for (i = 0; i < enspan.size(); i++)
  {
    it = key.find(enspan[i].cluster_id);
    if (it->second >= min_reads)   ////only consider cluster which contains more than n pairs
    {
      it2 = key_cl.find(enspan[i].cluster_id);//// key_cl contains  info that
      if (it2 == key_cl.end())
      {
        k++;
        enspan[i].cluster = k;
        key_cl[enspan[i].cluster_id] = k;
      }
      else
      {
        enspan[i].cluster = it2->second;
      }
      tmp_enspan.push_back(enspan[i]);
    }
  }
  enspan = tmp_enspan;
  return k;
}

/**
 * @brief write the programn parameters
 * @param inp_file
 * @param out_file
 * @param build
 * @param w
 * @param qual
 */
void write_enspan_params(string inp_file, string out_file, string build, double w, long qual)
{
  ofstream out;
  // Append parameters to parameter file
  out.open((out_file + "_params.txt").c_str());
  out << "ENSPAN" << std::endl;
  out << "inp_file\t" << inp_file << std::endl;
  out << "out_file\t" << out_file << std::endl;
  out << "qual\t" << qual << std::endl;
  out << "w\t" << w << std::endl;
  out << "build\t" << build << std::endl;
  out.close();
}

void write_enspan_out(string out_file, vector<cluster_info> &cluster, bool filter)
{
  bool condition_all, condition_filter;
  string hotspot, cosmic, hotspot_pair_match, cosmic_pair_match;
  sort(cluster.begin(), cluster.end(), cmp_cluster);
  ofstream out, out_filter;
  
  if (!filter)
  {
    out.open((out_file + "_fusion_all.txt").c_str());
    out << "Fusion_Type\t";
    out << "BreakPoint1\t";
    out << "BreakPoint2\t";
    out << "Gene1\tBreakPoint_Info_Pair1\t";
    out << "Gene2\tBreakPoint_Info_Pair2\t";
    out << "N_DRP\tN_SR\t";
    out << "BreakPoint1_Depth\tBreakPoint2_Depth\t";
    out << "BreakPoint1_AF\tBreakPoint2_AF\t";
    out << "BP1_Neighbour_Seq\tBP2_Neighbour_Seq\n";
  }
  
  out_filter.open((out_file + "_fusion.txt").c_str());
  out_filter << "Fusion_Type\t";
  out_filter << "BreakPoint1\t";
  out_filter << "BreakPoint2\t";
  out_filter << "Gene1\tBreakPoint_Info_Pair1\t";
  out_filter << "Gene2\tBreakPoint_Info_Pair2\t";
  out_filter << "N_DRP\tN_SR\t";
  out_filter << "BreakPoint1_Depth\tBreakPoint2_Depth\t";
  out_filter << "BreakPoint1_AF\tBreakPoint2_AF\t";
  out_filter << "BP1_Neighbour_Seq\tBP2_Neighbour_Seq\n";
  
  for (auto &i : cluster)
  {
    condition_filter =
      (i.n_split_read > 0 && i.p1_exact_pos != -1 && i.p2_exact_pos != -1) &&
      (!(i.p1_behalf_gene == "intergenic" && i.p2_behalf_gene == "intergenic") &&
       i.p1_behalf_gene != i.p2_behalf_gene) && (!i.is_rpt);
    condition_all = i.n_split_read > 0 && i.p1_exact_pos != -1 && i.p2_exact_pos != -1;
    
    hotspot = i.hotspot ? "1" : "0";
    cosmic = i.cosmic ? "1" : "0";
    hotspot_pair_match = i.sino_pair_match ? "1" : "0";
    cosmic_pair_match = i.cosmic_pair_match ? "1" : "0";
    
    if (condition_filter)
    {
      out_filter << i.fusion_type << "\t";
      out_filter << i.p1_chr << ":" << i.p1_exact_pos << "\t";
      out_filter << i.p2_chr << ":" << i.p2_exact_pos << "\t";
      out_filter << i.p1_behalf_gene << "\t";
      out_filter << i.p1_strand << ":" << i.p1_exon_info << "\t";
      out_filter << i.p2_behalf_gene << "\t";
      out_filter << i.p2_strand << ":" << i.p2_exon_info << "\t";
      out_filter << i.n_discordant_pair << "\t" << i.n_split_read << "\t";
      out_filter << i.p1_bp_depth << "\t" << i.p2_bp_depth << "\t";
      out_filter << i.p1_alle_freq << "\t" << i.p2_alle_freq << "\t";
      out_filter << i.p1_rpt << "\t" << i.p2_rpt << "\n";
      
    }
    
    if (!filter && condition_all)
    {
      out << i.fusion_type << "\t";
      out << i.p1_chr << ":" << i.p1_exact_pos << "\t";
      out << i.p2_chr << ":" << i.p2_exact_pos << "\t";
      out << i.p1_behalf_gene << "\t";
      out << i.p1_strand << ":" << i.p1_exon_info << "\t";
      out << i.p2_behalf_gene << "\t";
      out << i.p2_strand << ":" << i.p2_exon_info << "\t";
      out << i.n_discordant_pair << "\t" << i.n_split_read << "\t";
      out << i.p1_bp_depth << "\t" << i.p2_bp_depth << "\t";
      out << i.p1_alle_freq << "\t" << i.p2_alle_freq << "\t";
      out << i.p1_rpt << "\t" << i.p2_rpt << "\n";
    }
  }
  if (!filter)
    out.close();
  out_filter.close();
}

/**
 * @brief remove those clusters than are far away from others
 * @param enspans
 * @param w distance cutoff
 * @param rkey
 */
void remove_isolated_pairs(vector<discordant_pair> &enspans, double w)
{
  vector<discordant_pair> enspan_vec_tmp;
  sort(enspans.begin(), enspans.end(), cmp_p1_enspan_pairs);
  mask_pairs_chr_pos(enspans, w);
  if (!enspans.empty())
  {
    sort(enspans.begin(), enspans.end(), cmp_p2_enspan_pairs);
    mask_pairs_chr_pos(enspans, w);
    if (!enspans.empty())
    {
      sort(enspans.begin(), enspans.end(), cmp_p1_enspan_pairs);
    }
  }
}

void add_enspan_point_id(vector<discordant_pair> &enspan_vec)
{
  for (size_t i = 0; i < enspan_vec.size(); ++i)
  {
    enspan_vec[i].id = "pair_No_" + to_string(i);
  }
}

/**
 * @brief agglomerative hierarchical clustering method for enspan vector ,add cluster id for each of them
 * @param enspan_vec
 * @param distance_threshold
 * @param distance_type
 * @param out_matrix_file
 * @param min_reads_per_cluster
 * @return
 */
int
find_cluster_pairs_enspan_ahc(vector<discordant_pair> &enspan_vec, double distance_threshold, int distance_type,
                              int min_reads_per_cluster)
{
  vector<point> points;
  cluster_struct main_cluster;
  clock_t build_array_start, build_array_end, cluster_start, cluster_end;
  build_array_start = clock();
  build_pair_array(enspan_vec, points);//抽出坐标
  build_array_end = clock();
  cout << "the build cluster pair array step  costs time: " <<
       (build_array_end - build_array_start) / double(CLOCKS_PER_SEC) << " seconds " << std::endl;
  cluster_start = clock();
  init_cluster(main_cluster, distance_threshold, points, distance_type);
  cluster_end = clock();
  
  cout << "the init a cluster step  costs time: " <<
       (cluster_end - cluster_start) / double(CLOCKS_PER_SEC) << " seconds " << std::endl;
  add_cluster_id_for_enspan_vec(main_cluster, enspan_vec, min_reads_per_cluster);
  int root_cluster_num = print_root_nodes(main_cluster);
  return root_cluster_num;
}


void add_cluster_id_for_enspan_vec(cluster_struct &main_cluster,
                                   vector<discordant_pair> &enspan, int min_reads_per_cluster)
{
  vector<discordant_pair> tmp_enspan_vec;
  string cluster_id;
  int k = 0;
  for (int i = 0; i < main_cluster.num_nodes; ++i)
  {
    node temp_node = main_cluster.nodes[i];
    if (main_cluster.nodes[i].is_root && main_cluster.nodes[i].num_points >= min_reads_per_cluster)
    {
      cluster_id = "cluster_No_" + std::to_string(k);
      for (int j = 0; j < main_cluster.nodes[i].num_points; ++j)
      {
        int reserved_point_index = main_cluster.nodes[i].points[j];
        enspan[reserved_point_index].cluster_id = cluster_id;
        enspan[reserved_point_index].cluster = k;
        tmp_enspan_vec.push_back(enspan[reserved_point_index]);
      }
      k++;
    }
  }
  enspan = tmp_enspan_vec;
  tmp_enspan_vec.clear();
}

/**
 * @brief scan potential spanning pairs
 * @param inp_file
 * @param enspan
 * @param build
 * @param qual
 * @param w
 */
void
scan_discordant_pairs(const string &inp_file, const string &build, long qual, double w,
                      map<string, vector<discordant_pair>> &enspan_map, string nib_dir)
{
  //io-streams
  ifstream in;
  
  //variables
  //bool is_eof;
  string ref_name, tmp;
  vector<discordant_pair> enspan;
  
  //reference names
  map<string, string> key_ref;
  map<string, string>::iterator it_ref;
  
  //variables for misplaced reads
  map<string, bam_reduced> readname_2_alignment;
  map<string, bam_reduced>::iterator it_mpr;
  bam_reduced tmp_mpr;
  
  //bam variables
  samfile_t *fp;
  bam1_t *b = bam_init1();
  bam_map m1;
  
  //output
  discordant_pair tmp_enspan;
  //open bamfile
  fp = samopen(inp_file.c_str(), "rb", 0);
  if (fp == NULL)
  {
    cerr << "Error: cannot open " << inp_file << ".\n";
    exit(1);
  }
  
  //read reference list
  in.open((nib_dir + "/ref_names.txt").c_str());
  if (!in.is_open())
  {
    cerr << "Error: cannot open reference names file.\n";
    exit(1);
  }
  
  while (in >> tmp)
    key_ref[tmp] = tmp;
  in.close();
  
  bam_header_t *header = fp->header;
  
  cout << "Scanning discordant read pairs ...\n";
  
  while (samread(fp, b) > 0)
  {
    bam1_core_t c = b->core;
    /// if this read is mapped and paired, but not proper paired
    /// also not PCR duplicated and not secondary alignment
    if (c.qual >= qual && !(c.flag & BAM_FDUP) && !(c.flag & BAM_FSECONDARY) &&
        (c.flag & BAM_FPAIRED) && !(c.flag & BAM_FPROPER_PAIR))
    {
      read_bam_reduced_record(b, m1, header);
      //search for mate
      it_mpr = readname_2_alignment.find(m1.qname);
      /// find the mate in the map
      if (it_mpr != readname_2_alignment.end())
      {
        if (it_mpr->second.rname != m1.rname || abs(m1.pos - it_mpr->second.pos) >= w)
        {
          //index pairs
          uint32_t p1_chr_pos = combine_genome_chr_pos(header, b->core.tid, b->core.pos);
          uint32_t p2_chr_pos = combine_genome_chr_pos(header, b->core.mtid, b->core.mpos);
          //output misplaced mate
          if (p1_chr_pos <= p2_chr_pos)
          {
            tmp_enspan.qname = m1.qname;
            tmp_enspan.p1_flag = m1.flag;
            tmp_enspan.p1_chr = m1.rname;
            tmp_enspan.p1_pos = (uint32_t) m1.pos;
            tmp_enspan.p1_mapq = m1.mapq;
            tmp_enspan.p1_chr_pos = p1_chr_pos;
            tmp_enspan.p2_chr_pos = p2_chr_pos;
            tmp_enspan.p2_flag = it_mpr->second.flag;
            tmp_enspan.p2_chr = it_mpr->second.rname;
            tmp_enspan.p2_pos = (uint32_t) it_mpr->second.pos;
            tmp_enspan.p2_mapq = it_mpr->second.mapq;
            tmp_enspan.is_isolated = 0;
          }
          else
          {
            tmp_enspan.qname = m1.qname;
            tmp_enspan.p2_flag = m1.flag;
            tmp_enspan.p2_chr = m1.rname;
            tmp_enspan.p2_pos = (uint32_t) m1.pos;
            tmp_enspan.p2_mapq = m1.mapq;
            tmp_enspan.p1_chr_pos = p2_chr_pos;
            
            tmp_enspan.p2_chr_pos = p1_chr_pos;
            
            tmp_enspan.p1_flag = it_mpr->second.flag;
            tmp_enspan.p1_chr = it_mpr->second.rname;
            tmp_enspan.p1_pos = (uint32_t) it_mpr->second.pos;
            tmp_enspan.p1_mapq = it_mpr->second.mapq;
            tmp_enspan.is_isolated = 0;
          }
          //strand assignment for the pairs
          if (tmp_enspan.p1_flag & BAM_FREVERSE)
          {
            tmp_enspan.p1_strand = '-';
          }
          
          else
            tmp_enspan.p1_strand = '+';
          
          if (tmp_enspan.p2_flag & BAM_FREVERSE)
            tmp_enspan.p2_strand = '-';
          else
            tmp_enspan.p2_strand = '+';
          
          enspan.push_back(tmp_enspan);
        }
        readname_2_alignment.erase(it_mpr);
      }
        /// the first mate
      else
      {
        //add read to buffer
        tmp_mpr.qname = m1.qname;
        tmp_mpr.flag = m1.flag;
        tmp_mpr.rname = m1.rname;
        tmp_mpr.pos = m1.pos;
        tmp_mpr.mapq = m1.mapq;
        readname_2_alignment[m1.qname] = tmp_mpr;
      }
    }
  }
  cout << "Scanning discordant read pairs done.\n";
  samclose(fp);
  
  set<string> chr_pair_set;
  get_chr_pair_set(enspan, chr_pair_set);
  vector<discordant_pair> tmp_enspan_vec;
  
  for (auto set_it = chr_pair_set.begin(); set_it != chr_pair_set.end(); ++set_it)
  {
    enspan_map[*set_it] = tmp_enspan_vec;
  }
  
  for (auto &i : enspan)
  {
    enspan_map[i.p1_chr + "_" + i.p2_chr].push_back(i);
  }
  
  
}

/**
 *
 * @param txpts
 * @param chr1
 * @param chr2
 * @param p1_pos
 * @param p2_pos
 * @param exon_infos
 * @param p1_part
 * @param p2_part
 */
void
add_exon_anno(vector<RefSeqTranscript> &txpts, string chr1, string chr2, long &p1_pos, long &p2_pos,
              vector<string> &exon_infos, const string &p1_part, const string &p2_part)
{
  int record_count_p1;
  int record_count_p2;
  record_count_p1 = 0;
  record_count_p2 = 0;
  string p1_gene, p2_gene, p1_txpt, p2_txpt, p1_exon_info, p2_exon_info, p1_strand, p2_strand;
  p1_gene = p2_gene = p1_txpt = p2_txpt = p1_exon_info = p2_exon_info = p1_strand = p2_strand = "";
  vector<int> p1_exon_no, p2_exon_no;
  RefSeqTranscript tmp_txpt_p1, tmp_txpt_p2;
  vector<RefSeqTranscript> tmp_txpts_vec_p1, tmp_txpts_vec_p2;
  string p1_gene_vec_str, p2_gene_vec_str;
  set<string> p1_gene_vec, p2_gene_vec;
  set<string>::iterator genes_it;
  string p1_gene_part = "";
  string p2_gene_part = "";
  string n_terminal_gene = "";
  string c_terminal_gene = "";
  string fusion_5_3_exon_pair;
  if (p1_pos != -1)
  {
    tmp_txpts_vec_p1.clear();
    for (size_t i = 0; i < txpts.size(); ++i)
    {
      if (record_count_p1 == 0 &&
          (chr1 == txpts[i].chrom && p1_pos >= txpts[i].txStart &&
           p1_pos <= txpts[i].txEnd))
      {
        tmp_txpts_vec_p1.push_back(txpts[i]);
      }
    }
    // if p1 can not find annotate
    if (tmp_txpts_vec_p1.size() == 0)
    {
      p1_exon_info = ".";
      p1_gene = "intergenic";
      p1_strand = ".";
    }
    else
    {
      find_genes_from_txpts(tmp_txpts_vec_p1, p1_gene_vec);
      find_the_longest_cds_txpt(tmp_txpts_vec_p1, tmp_txpt_p1);
      p1_gene = tmp_txpt_p1.geneName;
      p1_txpt = tmp_txpt_p1.transcriptID;
      p1_strand = tmp_txpt_p1.strand;
      add_exon_num_anno(tmp_txpt_p1, p1_pos, chr1, p1_exon_no);
      p1_exon_info = p1_txpt + ":" + to_string(p1_exon_no[0]) + "-" + to_string(p1_exon_no[1]);
    }
  }
  else
  {
    p1_exon_info = ".";
    p1_gene = ".";
    p1_strand = ".";
    
  }
  
  if (p2_pos != -1)
  {
    tmp_txpts_vec_p2.clear();
    for (size_t i = 0; i < txpts.size(); ++i)
    {
      if (record_count_p2 == 0 &&
          (chr2 == txpts[i].chrom && p2_pos >= txpts[i].txStart &&
           p2_pos <= txpts[i].txEnd))
      {
        tmp_txpts_vec_p2.push_back(txpts[i]);
      }
    }
    //if p2 can not find annotate
    if (tmp_txpts_vec_p2.size() == 0)
    {
      p2_exon_info = ".";
      p2_gene = "intergenic";
      p2_strand = ".";
    }
    else
    {
      find_genes_from_txpts(tmp_txpts_vec_p2, p2_gene_vec);
      find_the_longest_cds_txpt(tmp_txpts_vec_p2, tmp_txpt_p2);
      p2_gene = tmp_txpt_p2.geneName;
      p2_txpt = tmp_txpt_p2.transcriptID;
      p2_strand = tmp_txpt_p2.strand;
      add_exon_num_anno(tmp_txpt_p2, p2_pos, chr2, p2_exon_no);
      p2_exon_info = p2_txpt + ":" + to_string(p2_exon_no[0]) + "-" + to_string(p2_exon_no[1]);
    }
  }
  else
  {
    p2_exon_info = ".";
    p2_gene = ".";
    p2_strand = ".";
  }
  
  p1_gene_vec_str = "";
  p2_gene_vec_str = "";
  if (p1_gene_vec.size() > 0)
  {
    for (genes_it = p1_gene_vec.begin(); genes_it != p1_gene_vec.end(); ++genes_it)
    {
      if (p1_gene_vec_str == "")
        p1_gene_vec_str = (string) genes_it->data();
      else
        p1_gene_vec_str = p1_gene_vec_str + ";" + (string) genes_it->data();
    }
  }
  else
  {
    p1_gene_vec_str = ".";
  }
  
  if (p2_gene_vec.size() > 0)
  {
    for (genes_it = p2_gene_vec.begin(); genes_it != p2_gene_vec.end(); ++genes_it)
    {
      if (p2_gene_vec_str == "")
        p2_gene_vec_str = (string) genes_it->data();
      else
        p2_gene_vec_str = p2_gene_vec_str + ";" + (string) genes_it->data();
    }
  }
  else
  {
    p2_gene_vec_str = ".";
  }
  
  
  //TODO: add gene part info
  
  string p1_bp_exon, p2_bp_exon;
  if (p1_strand == ".")
  {
    p1_gene_part = ".";
    p1_bp_exon = ".";
  }
  else
  {
    if ((p1_strand == "+" && p1_part == "left") || (p1_strand == "-" && p1_part == "right"))
    {
      p1_gene_part = "upstream";
      p1_bp_exon = to_string(p1_exon_no[0]);
    }
    
    if ((p1_strand == "+" && p1_part == "right") || (p1_strand == "-" && p1_part == "left"))
    {
      p1_gene_part = "downstream";
      p1_bp_exon = to_string(p1_exon_no[1]);
    }
  }
  
  if (p2_strand == ".")
  {
    p2_gene_part = ".";
    p2_bp_exon = ".";
  }
  else
  {
    if ((p2_strand == "+" && p2_part == "left") || (p2_strand == "-" && p2_part == "right"))
    {
      p2_gene_part = "upstream";
      p2_bp_exon = to_string(p2_exon_no[0]);
    }
    
    if ((p2_strand == "+" && p2_part == "right") || (p2_strand == "-" && p2_part == "left"))
    {
      p2_gene_part = "downstream";
      p2_bp_exon = to_string(p2_exon_no[1]);
    }
  }
  
  string fusion_5exon_pair, fusion_3exon_pair, fusion_pair;
  string up_gene = "", down_gene = "";
  fusion_5_3_exon_pair = "";
  fusion_5exon_pair = "";
  fusion_3exon_pair = "";
  fusion_pair = "";
  if ((p1_gene != "intergenic" && p2_gene != "intergenic") && (p1_gene_part != p2_gene_part))
  {
    if (p1_gene_part == "upstream")
    {
      up_gene = p1_gene;
      down_gene = p2_gene;
      fusion_pair = p1_gene + "," + p2_gene;
    }
    else
    {
      up_gene = p2_gene;
      down_gene = p1_gene;
      fusion_pair = p2_gene + "," + p1_gene;
    }
  }
  else
  {
    up_gene = ".";
    down_gene = ".";
    fusion_pair = ".,.";
  }
  exon_infos.push_back(p1_gene);
  exon_infos.push_back(p2_gene);
  exon_infos.push_back(p1_exon_info);
  exon_infos.push_back(p2_exon_info);
  exon_infos.push_back(p1_strand);
  exon_infos.push_back(p2_strand);
  exon_infos.push_back(p1_gene_vec_str);
  exon_infos.push_back(p2_gene_vec_str);
  exon_infos.push_back(p1_gene_part);
  exon_infos.push_back(p2_gene_part);
  exon_infos.push_back(p1_bp_exon);
  exon_infos.push_back(p2_bp_exon);
  exon_infos.push_back(up_gene);
  exon_infos.push_back(down_gene);
  exon_infos.push_back(fusion_pair);
  
  
}

/**
 * @brief add the exon number information based on refseq gene annotation
 * @param txpt
 * @param pos
 * @param chr
 * @param exon_nums
 */
void add_exon_num_anno(RefSeqTranscript &txpt, long &pos, string &chr, vector<int> &exon_nums)
{
  int start_exon_no, end_exon_no, temp_index;
  start_exon_no = end_exon_no = 0;
  for (int i = 0; i < txpt.codingExonParts.size() - 1; ++i)
  {
    if (pos >= txpt.codingExonParts[i] && pos <= txpt.codingExonParts[i + 1])
    {
      temp_index = i / 2 + 1;
      if (txpt.strand == "+")
      {
        if (i % 2 == 1)//2n+1
        {
          start_exon_no = temp_index;
          end_exon_no = temp_index + 1;
        }
        else
        {
          start_exon_no = temp_index;
          end_exon_no = temp_index;
        }
      }
      if (txpt.strand == "-")
      {
        if (i % 2 == 1)//奇数
        {
          start_exon_no = txpt.codingExonCout + 1 - (temp_index + 1);
          end_exon_no = txpt.codingExonCout + 1 - temp_index;
        }
        else
        {
          start_exon_no = txpt.codingExonCout + 1 - (temp_index + 1);
          end_exon_no = txpt.codingExonCout + 1 - (temp_index + 1);
        }
      }
      break;
    }
  }
  exon_nums.push_back(start_exon_no);
  exon_nums.push_back(end_exon_no);
}

void build_pair_array(vector<discordant_pair> &enspan_vec, vector<point> &points)//传递指向points数组的指针
{
  for (size_t i = 0; i < enspan_vec.size(); ++i)
  {
    point p;
    p.pos.x = enspan_vec[i].p1_chr_pos;
    p.pos.y = enspan_vec[i].p2_chr_pos;
    p.label = "x=" + to_string(enspan_vec[i].p1_chr_pos) + ",y=" +
              to_string(enspan_vec[i].p2_chr_pos);
    points.push_back(p);
  }
}

/**
 * @brief check the inter-distance between 2 clusters
 * @param enspan
 * @param distance
 */
void mask_pairs_chr_pos(vector<discordant_pair> &enspan, long distance)
{
  long Lx, Ly;
  long lr, ll, np;
  vector<discordant_pair> enspans_tmp;
  //long i;
  //int n = distance;
  np = enspan.size();
  if (np <= 2)
  {
    enspan.clear();
    return;
  }
  
  if (np >= 3)
  {
    //first read pair
    Lx = abs((int32_t) (enspan[1].p1_chr_pos - enspan[2].p1_chr_pos));
    Ly = abs((int32_t) (enspan[1].p2_chr_pos - enspan[2].p2_chr_pos));
    if (!(Lx > distance || Ly > distance))
    {
      enspans_tmp.push_back(enspan[1]);
    }
    else
    {
      enspan[1].is_isolated = 1;
    }
    //last read pair
    Lx = abs((int32_t) (enspan[np - 1].p1_chr_pos - enspan[np - 2].p1_chr_pos));
    Ly = abs((int32_t) (enspan[np - 1].p2_chr_pos - enspan[np - 2].p2_chr_pos));
    if (Lx > distance || Ly > distance)
      enspan[np - 1].is_isolated = 1;
    for (int i = 1; i < np - 1; i++)
    {
      ll = abs((int32_t) (enspan[i - 1].p1_chr_pos - enspan[i].p1_chr_pos));
      lr = abs((int32_t) (enspan[i + 1].p1_chr_pos - enspan[i].p1_chr_pos));
      if (ll < lr)
        Lx = ll;
      else
        Lx = lr;
      ll = abs((int32_t) (enspan[i - 1].p2_chr_pos - enspan[i].p2_chr_pos));
      lr = abs((int32_t) (enspan[i + 1].p2_chr_pos - enspan[i].p2_chr_pos));
      if (ll < lr)
        Ly = ll;
      else
        Ly = lr;
      
      if (!(Lx > distance || Ly > distance))
      {
        enspans_tmp.push_back(enspan[i]);
      }
      else
      {
        enspan[i].is_isolated = 1;
      }
//    if (Lx > distance || Ly > distance)
//      enspan[i].is_isolated = 1;
    }
    
    enspan = enspans_tmp;
    return;
  }
  
  
}

void get_chr_pair_set(vector<discordant_pair> &enspan, set<string> &chr_pair_set)
{
  
  for (int i = 0; i < enspan.size(); ++i)
  {
    chr_pair_set.insert(enspan[i].p1_chr + "_" + enspan[i].p2_chr);
  }
}

string determine_fusion_type_from_drp(cluster_info &cluster)
{
  string fusion_type = "";
  auto iterator1 = cluster.drp_type_set.find("diff_chr");
  auto iterator2 = cluster.drp_type_set.find("same_chr_with_same_orientation");
  auto iterator3 = cluster.drp_type_set.find("same_chr_with_absolute_reverse");
  auto iterator4 = cluster.drp_type_set.find("same_chr_with_default_orientation");
  
  if (iterator1 != cluster.drp_type_set.end())
    fusion_type = "Translocation";
  if (iterator2 != cluster.drp_type_set.end())
    fusion_type = "Inversion";
  if (iterator3 != cluster.drp_type_set.end())
    fusion_type = "Duplication";
  if (iterator4 != cluster.drp_type_set.end())
    fusion_type = "Deletion";
  if (fusion_type == "")
    fusion_type = "Unknown";
  return fusion_type;
}

void get_mean_insert_size(string input_bam, vector<double> &insert)
{
  vector<int> insert_size;
  long insert_size_total = 0;
  long insert_size_sd_total = 0;
  uint32_t filter_flag = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
  samFile *fp;
  /// read the bam file
  fp = samtools_sam_open(input_bam.c_str());
  if (fp == NULL)
  {
    std::cerr << "Error: can not open bam-file: " << input_bam << std::endl;
    exit(1);
  }
  bam1_t *b;
  bam1_core_t *c;
  b = bam_init1();
  c = &b->core;
  int ret;
  bam_hdr_t *header = sam_hdr_read(fp);
  while ((ret = sam_read1(fp, header, b)) >= 0)
  {
    
    if (((c->flag & BAM_FPAIRED) && (c->flag & BAM_FPROPER_PAIR)) && (!(c->flag & filter_flag)))
      insert_size.push_back(abs(c->isize));
  }
  
  vector<int>::iterator iter;
  for (iter = insert_size.begin(); iter != insert_size.end(); iter++)
  {
    insert_size_total += *iter;
  }
  double insert_size_mean = (double) insert_size_total / (double) insert_size.size();
  for (iter = insert_size.begin(); iter != insert_size.end(); iter++)
  {
    insert_size_sd_total += (((double) (*iter)) - insert_size_mean) * (((double) (*iter)) - insert_size_mean);
  }
  double insert_size_sd = sqrt(insert_size_sd_total / (double) insert_size.size());
  bam_destroy1(b);
  bam_hdr_destroy(header);
  sam_close(fp);
  insert.clear();
  insert.push_back(insert_size_mean);
  insert.push_back(insert_size_sd);
  
}