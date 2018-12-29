#pragma once

#include <stdint.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <regex>

#include "nibtools.h"
#include "installdir.h"
#include "util_bam.h"
using namespace std;

class RefSeqTranscript
{
public:
  RefSeqTranscript();
  
  RefSeqTranscript(std::string refseq_line);
  
  ~RefSeqTranscript();
  
  ///  split string into int vector
  std::vector<uint32_t> splitStringToInt(const std::string &s, const std::string &delim);
  
  /// remove UTR from exons
  void removeUTR();
  /// find the longest txpt for gene
 
  
  
  void printOut(std::ostream &os);
  /// start position is 0-based inclusive
  /// end position is 0-based exclusive, ie, [start,end)
  std::string bin; /// bin field
  std::string transcriptID;          /// "Name of gene (usually transcript_id from GTF)"
  std::string chrom;        /// "Chromosome name"
  std::string strand;      /// "+ or - for strand"
  uint32_t txStart;        /// "Transcription start position"
  uint32_t txEnd;          /// "Transcription end position"
  uint32_t cdsStart;        /// "Coding region start"
  uint32_t cdsEnd;          /// "Coding region end"
  uint32_t exonCount;      /// "Number of exons"
  std::vector<uint32_t> exonStarts; /// "Exon start positions", including 5' UTR
  std::vector<uint32_t> exonEnds;   /// "Exon end positions", including 3' UTR
  int score;              /// "Score"
  std::string geneName;        /// "Alternate name (e.g. gene_id from GTF)"
  std::string cdsStartStat;  /// "enum('none','unk','incmpl','cmpl')"
  std::string cdsEndStat;    /// "enum('none','unk','incmpl','cmpl')"
  std::string exonFrames;  /// "Exon frame offsets {0,1,2}"
  
  std::vector<uint32_t> codingExonStarts; /// "coding Exon start positions, UTR will be excluded"
  std::vector<uint32_t> codingExonEnds;   /// "coding Exon end positions, UTR will be excluded"
  int codingExonCout;
  std::vector<uint32_t> codingExonParts;
  /// from 3' UTR to 5' UTR
  uint32_t transcriptLength;
  
  /// transcribed cDNA length
  uint32_t cDNALength;

private:
  
};

void
readRefSeqTranscript(std::string filename, std::vector<RefSeqTranscript> &refSeqTranscripts);

void
find_the_longest_txpts_per_gene(vector<RefSeqTranscript> &txpts, vector<RefSeqTranscript> &longest_txpts);

void find_the_longest_cds_txpt(vector<RefSeqTranscript> &txpts,RefSeqTranscript& longestCdsTxpt);
void find_genes_from_txpts(vector<RefSeqTranscript> &txpts_vec, set<string> &genes_vec);

void add_cds_parts(vector<RefSeqTranscript> &txpts);
