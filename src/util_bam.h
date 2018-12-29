#pragma once

#include <string>
#include <sam.h>
#include <htslib/faidx.h>
#include <set>
#include <sstream>
#include <htslib/tbx.h>
#include <vector>


using namespace std;
#define MAXCIGAR 1024



typedef struct
{
  std::string qname;
  long flag;
  std::string rname;
  long pos;
  long mapq;
  int cigar_pos[MAXCIGAR];
  char cigar_type[MAXCIGAR];
  int n_cigar;
  std::string mrnm;
  long mpos;
  long isize;
  std::vector<int> seq;
  std::string qual;
  long size;
  std::vector<int> mask_subst;
  int mask_size;
  bool is_indel;
  bool is_MID;
} bam_map;


typedef struct
{
  long index;
  std::string qname;
  long flag;
  std::string rname;
  long pos;
  long mapq;
  std::string mrnm;
  long mpos;
  long isize;
} bam_reduced;




void read_bam_reduced_record( bam1_t *b, bam_map &m, bam_header_t *header);
uint32_t combine_genome_chr_pos(bam_header_t *header, int chromID, int32_t position);
std::string get_right_neighbor_sequence_nib(std::string chrom, int32_t pos_1based, int length, string nib);
std::string get_left_neighbor_sequence_nib(std::string chrom, int32_t pos_1based, int length, string nib);
std::string get_sequence_nib(std::string chrom, int32_t start_1based, int32_t end_1based, string nib);
std::string chromID2ChrName(int refID);



