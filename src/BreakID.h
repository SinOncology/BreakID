#pragma once

#include <vector>
#include <fstream>
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
#include <getopt.h>
#include <cmath>

#include "util_bam.h"
#include "util_bed.h"
#include "util_cluster.h"
#include "installdir.h"
#include "RefSeqTranscript.h"
#include "BamAlignment.h"

string BreakID_help = " Usage: \n \t BreakID -i input.bam -o prefix -n nib_folder <options> \n\n \
     DESCRIPTION\n \
     \t -h -? -help \t help\n \
     \t -i*        \t input bam-file\n \
     \t -o*        \t output file (prefix only)\n \
     \t -n*        \t folder name to nib files\n \
     \t -q         \t encompassing reads quality thresholds  [20]\n\
     \t -t         \t distance relative to (sqrt(2)*(insert size mean +3* insert size sd))  [2]\n \
     \t -fast      \t use the fast cluster strategy [default no] \n \
     \t -all       \t no filter enspan out [default is filter]  \n ";


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
  set <string> drp_type_set;
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


/// prototypes
void scan_discordant_pairs(const string &inp_file, const string &build, long qual, double w,
                           map <string, vector<discordant_pair>> &enspan_map, string nib_dir);

void annotate_cluster_for_sa_tag(vector <cluster_info> &clusters, string nib_dir);

void get_mean_insert_size(string input_bam, vector<double> &insert);

void get_chr_pair_set(vector <discordant_pair> &enspan, set <string> &chr_pair_set);

void mask_pairs_chr_pos(vector <discordant_pair> &enspan, long distance);

void build_pair_array(vector <discordant_pair> &, vector <point> &);

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
find_cluster_pairs_enspan_ahc(vector <discordant_pair> &enspan, double distance_threshold, int distance_type,
                              int min_reads_per_cluster);

void
add_cluster_id_for_enspan_vec(cluster_struct &main_cluster, vector <discordant_pair> &enspan,
                              int min_reads_per_cluster);

void add_enspan_point_id(vector <discordant_pair> &enspan_vec);

void remove_isolated_pairs(vector <discordant_pair> &enspans, double w);

void write_enspan_out(string out_file, vector <cluster_info> &cluster, bool filter);

void write_enspan_params(string inp_file, string out_file, string build, double w, long qual);

void
add_exon_anno(vector <RefSeqTranscript> &txpts, string chr1, string chr2, long &p1_pos, long &p2_pos,
              vector <string> &exon_infos, const string &p1_part, const string &p2_part);

void add_exon_num_anno(RefSeqTranscript &txpt, long &pos, string &chr, vector<int> &exon_nums);

int find_cluster_pairs_enspan_fast(vector <discordant_pair> &enspan, double w, int min_reads);

void
findClusterBreakPointInfoSaTag(string bam_file, vector <discordant_pair> &enspan, double w,
                               vector <cluster_info> &cluster_vec, vector<bam1_t *> &split_reads, string nib_dir);

void find_sa_reads(samfile_t *fp, const string region_chr, uint32_t region_start, uint32_t region_end,
                   map <string, vector<split_align_pair>> &encompassing_map, bam_index_t *idx_bam);

void findEncompassingReadsAndBreakPointInfo(map<long, cluster_info> &cluster_map, const string bam_file, const int w,
                                            vector<bam1_t *> &split_reads);

void find_bp_pair(map <string, vector<split_align_pair>> &p1_encompass_map,
                  map <string, vector<split_align_pair>> &p2_encompass_map, breakpoint_pair &bp_pair,
                  const string &p1_chr, const string &p2_chr, vector<bam1_t *> &split_reads, int bp_pos_error);

string determine_fusion_type_from_drp(cluster_info &cluster);
