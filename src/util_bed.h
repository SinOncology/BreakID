#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <algorithm>
#include <sam.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>


using namespace std;




double cal_single_base_depth(string chrom, uint64_t pos, samfile_t *fp, bam_index_t *idx_bam);

double cal_mean_depth_oc(std::string chrom, uint64_t start, uint64_t end, samfile_t *fp, bam_index_t *idx_bam);

double cal_mean_depth(std::string chrom, uint64_t start, uint64_t end, samfile_t *fp, bam_index_t *idx_bam);

std::vector<std::string> split_string(const std::string &s,
                                      const std::string &delim);

int find_longest_repeat_substring(const std::string &s);

typedef struct
{
  std::string rawstring;
  uint16_t start_index;
  std::string sub_str;
  uint16_t length;
} repeat_str;