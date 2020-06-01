//#include "util.h"
#include <sam.h>
#include "util_bam.h"
#include "nibtools.h"


void read_bam_reduced_record( bam1_t *b, bam_map &m, bam_header_t *header)
{
 
    std::stringstream line;
  
    const bam1_core_t *c = &b->core;
    
    m.qname = (std::string) bam1_qname(b);
    m.flag = (long) c->flag;
    
    if (c->tid < 0)
      m.rname = "*";
    else
      m.rname = header->target_name[c->tid];
    
    m.pos = c->pos + 1;
    m.mapq = c->qual;
    //mrnm
    if (c->mtid < 0)
      m.mrnm = "*";
    else if (c->mtid == c->tid)
      m.mrnm = "=";
    else
    {
      if (header)
        m.mrnm = header->target_name[c->mtid];
      else
      {
        line.str("");
        line.clear();
        line << c->mtid;
        m.mrnm = line.str();
      }
    }
    
    m.mpos = c->mpos + 1;
    m.isize = c->isize;
    m.size = c->l_qseq;
  
  
}


/**
 * @brief given a 0-based position in a chromosome, return a 0-based genome-wide position
 * @param header bam file header
 * @param chromID 0-based chromsome ID
 * @param position 0-based position
 * @return 0-based genome-wide position
 */
uint32_t combine_genome_chr_pos(bam_header_t *header, int chromID, int32_t position)
{
  /// chr_pos is 0-based
  uint32_t chr_pos = 0;
  for (int i = 0; i < chromID; ++i)
  {
   chr_pos += header->target_len[i];
  }
  /// position is 0-based, chr_pos is also 0-based
  chr_pos  += position;
  return chr_pos;
}


/**
 * @brief given a position, get the bases (length) in the right neighbor (do not include the base in the given position)
 * @param chrom chrom name
 * @param pos_1based 1-based reference position, for example, from vcf format
 * @param length, the length of the wanted bases
 * @return
 */
std::string get_right_neighbor_sequence_nib(std::string chrom, int32_t pos_1based, int length, string nib_dir)
{
  std::string seq="";
  //nib access
  nib nibObj;
  std::string build="hg19";
 
  //open genome access (nib)
  std::string fn= nib_dir + "/" + build + "_" + chrom + ".nib";
  int stat=nibObj.open(fn);
  char base;
  for (int32_t i = pos_1based; i< pos_1based + length; i++)
  {
    stat = nibObj.getBase(&base, i);
    seq +=base;
  }
  nibObj.close();
  return seq;
}

/**
 * @brief given a positition, get the bases (length) in the left neighbor (do not include the base in the given position)
 * @param chrom chrom name
 * @param pos_1based 1-based reference position, for example, from vcf format
 * @param length, the length of the wanted bases
 * @return
 */
std::string get_left_neighbor_sequence_nib(std::string chrom, int32_t pos_1based, int length, string nib_dir)
{
  std::string seq="";
  //nib access
  nib nibObj;
  std::string build="hg19";
  //open genome access (nib)
  std::string fn= nib_dir + "/" + build + "_" + chrom + ".nib";
  int stat=nibObj.open(fn);
  char base;
  for (int32_t i = pos_1based -length; i< pos_1based ; i++)
  {
    stat = nibObj.getBase(&base, i-1);
    seq +=base;
  }
  nibObj.close();
  return seq;
}




/// refID is 0-based in bam file
std::string chromID2ChrName(int refID)
{
  std::string chromName = "";
  std::ostringstream os;
  if (refID == 23)
    chromName = "chrY";
  else if (refID == 22)
    chromName = "chrX";
  else if (refID >= 0 && refID < 22)
  {
    os << "chr" << refID + 1;
    chromName = os.str();
  }
  return (chromName);
}




