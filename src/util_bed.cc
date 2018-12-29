#include <htslib/sam.h>
#include "util_bed.h"
#include "BamAlignment.h"





/**
 * @brief calculate the mean depth of a given target region
 * @param chrom
 * @param start
 * @param end
 * @param fp point to bam file
 * @param idx_bam  pointer to bai file
 * @return
 */
double cal_mean_depth(string chrom, uint64_t start, uint64_t end, samfile_t *fp, bam_index_t *idx_bam)
{
  /// Parse region
  int tid_t;
  /// chromosome start and end position
  int beg_t, end_t;
  uint32_t n_bases = end - start + 1;
  long coverage = 0;
  double mean_coverage = 0;
  bam_parse_region(fp->header, chrom.c_str(), &tid_t, &beg_t, &end_t);
  
  bam_iter_t iter_bam;
  
  iter_bam = bam_iter_query(idx_bam, tid_t, start, end);
  
  bam1_t *b = bam_init1();
  
  
  uint32_t rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
  /// Iterate every alignment in this region
  while (bam_iter_read(fp->x.bam, iter_bam, b) >= 0)
  {
    if (!(b->core.flag & rflag_filter))
    {
      BamAlignment alignment(b);
      long bam_start = alignment.getAlignmentStart();
      long bam_end = alignment.getAlignmentEnd();
      
      // this alignment does not coverage start-end
      if (bam_end < start || bam_start > end)
        continue;
      
      if (bam_start <= start)
      {
        if (bam_end <= end)
          coverage += bam_end - start + 1;
        else
          coverage += end - start + 1;
      }
      else
      {
        if (bam_end <= end)
          coverage += bam_end - bam_start + 1;
        else
          coverage += end - bam_start + 1;
      }
      
    }
  }
  mean_coverage = (double) coverage / (double) n_bases;
  bam_destroy1(b);
  return mean_coverage;
}

/**
 * @brief calculate the mean depth for a given target region, using the OC tag
 * @param chrom
 * @param start
 * @param end
 * @param fp
 * @param idx_bam
 * @return
 */
double cal_mean_depth_oc(string chrom, uint64_t start, uint64_t end, samfile_t *fp, bam_index_t *idx_bam)
{
  /// Parse region
  int tid_t;
  /// chromosome start and end position
  int beg_t, end_t;
  uint32_t n_bases = end - start + 1;
  long coverage = 0;
  double mean_coverage = 0;
  long bam_start = 0;
  long bam_end = 0;
  bam_parse_region(fp->header, chrom.c_str(), &tid_t, &beg_t, &end_t);
  bam1_t *b = bam_init1();
  
  bam_iter_t iter_bam;
  samread(fp, b);
  int read_length = b->core.l_qseq;
  
  iter_bam = bam_iter_query(idx_bam, tid_t, start - read_length, end);
  
  
  uint32_t rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
  /// Iterate every alignment in this region
  while (bam_iter_read(fp->x.bam, iter_bam, b) >= 0)
  {
    if (!(b->core.flag & rflag_filter))
    {
      BamAlignment alignment(b);
      string oc_tag = alignment.originCigar();
      bam_start = alignment.getAlignmentStart();
      if (oc_tag != "")
      {
        CigarRoller tmp_cigar;
        tmp_cigar.Set(oc_tag.c_str());
        bam_end = (uint32_t) tmp_cigar.getAlignmentEnd(
          (uint32_t) alignment.getAlignmentStart());
      }
      else
        bam_end = alignment.getAlignmentEnd();
      
      // this alignment does not coverage start-end
      if (bam_end < start || bam_start > end)
        continue;
      
      if (bam_start <= start)
      {
        if (bam_end <= end)
          coverage += bam_end - start + 1;
        else
          coverage += end - start + 1;
      }
      else
      {
        if (bam_end <= end)
          coverage += bam_end - bam_start + 1;
        else
          coverage += end - bam_start + 1;
      }
    }
  }
  mean_coverage = (double) coverage / (double) n_bases;
  bam_destroy1(b);
  return mean_coverage;
}

/**
 * @brief calculate the depth for a given position
 * @param chromosome
 * @param pos
 * @param fp
 * @param idx_bam
 * @return
 */
double cal_single_base_depth(string chromosome, uint64_t pos, samfile_t *fp, bam_index_t *idx_bam)
{
  
  /// Parse region
  int tid_t;
  /// chromosome start and end position
  int beg_t, end_t;
  int depth = 0;
  
  
  bam_parse_region(fp->header, chromosome.c_str(), &tid_t, &beg_t, &end_t);
  if (pos < beg_t || pos >= end_t || pos > end_t)
  {
    std::cerr << "Error: region is beyond chromosome length" << std::endl;
    exit(1);
  }
  /// Create bam iterator
  bam_iter_t iter_bam;
  
  ///  bam_iter_query takes 0-based positions as parameters
  iter_bam = bam_iter_query(idx_bam, tid_t, pos - 1, pos);
  
  bam1_t *b = bam_init1();
  /// the filter flag for pileup, the same as the one in samtools mpileup
  
  bool condition;
  /// Iterate every alignment in this region
  while (bam_iter_read(fp->x.bam, iter_bam, b) >= 0)
  {
    condition = (b->core.qual > 0) && (!(b->core.flag & BAM_FDUP)) && (b->core.flag & BAM_FPAIRED);
    if (condition)
    {
      depth++;
    }
    
  } /// finish parsing all the reads in this region
  bam_destroy1(b);
  return depth;
}

std::vector<std::string> split_string(const std::string &s,
                                      const std::string &delim)
{
  bool keep_empty = false;
  std::vector<std::string> result;
  result.empty();
  if (delim.empty())
  {
    result.push_back(s);
    return result;
  }
  
  std::string::const_iterator substart = s.begin(), subend;
  while (true)
  {
    subend = search(substart, s.end(), delim.begin(), delim.end());
    std::string temp(substart, subend);
    if (keep_empty || !temp.empty())
    {
      result.push_back(temp);
    }
    if (subend == s.end())
    {
      break;
    }
    substart = subend + delim.size();
  }
  return result;
}

int find_longest_repeat_substring(const std::string &str)
{
  vector<repeat_str> repeat_vec;
  uint16_t len = str.length();
  uint16_t start_index = 0;
  repeat_str repeat_tmp;
  while (start_index < len)
  {
    uint16_t sub_len = 0;
    char begin_char = str.substr(start_index, 1).c_str()[0];
    bool same = true;
    while (same)
    {
      sub_len++;
      char current;
      current = str.substr(start_index + sub_len, 1).c_str()[0];
      same = (current == begin_char);
    }
    repeat_tmp.length = sub_len;
//    repeat_tmp.max_index = min_index+sub_len;
    repeat_tmp.start_index = start_index;
    repeat_tmp.rawstring = str;
    repeat_tmp.sub_str = str.substr(start_index, sub_len);
    repeat_vec.push_back(repeat_tmp);
    start_index = start_index + sub_len;
  }
  uint16_t max_len = 0;
  int max_repeat_index = 0;
  for (int i = 0; i < repeat_vec.size(); ++i)
  {
    if (repeat_vec[i].length > max_len)
    {
      max_len = repeat_vec[i].length;
      max_repeat_index = i;
    }
  }
  return repeat_vec[max_repeat_index].length;
}