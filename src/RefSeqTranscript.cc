
#include "RefSeqTranscript.h"



RefSeqTranscript::RefSeqTranscript()
{
  transcriptID = "";
  geneName = "";
  exonCount = 0;
  codingExonCout = 0;
  
}

/**
 * @brief parse a refGene line and assign the values
 * @param refseq_line
 */
RefSeqTranscript::RefSeqTranscript(std::string refseq_line)
{
  std::stringstream line;
  std::string tmp;
  line.str("");
  line.clear();
  line << refseq_line;
  getline(line, tmp, '\t');
  /// bin field
  bin = tmp;
  /// transcript ID
  getline(line, tmp, '\t');
  transcriptID = tmp;
  /// chrom name
  getline(line, tmp, '\t');
  chrom = tmp;
  /// strand
  getline(line, tmp, '\t');
  strand = tmp;
  /// transcript start
  getline(line, tmp, '\t');
  txStart = (uint32_t) std::atol(tmp.c_str());
  /// transcript end
  getline(line, tmp, '\t');
  txEnd = (uint32_t) std::atol(tmp.c_str());
  /// coding exon start
  getline(line, tmp, '\t');
  cdsStart = (uint32_t) std::atol(tmp.c_str());
  /// coding exon end
  getline(line, tmp, '\t');
  cdsEnd = (uint32_t) atol(tmp.c_str());
  /// number of exonx
  getline(line, tmp, '\t');
  exonCount = (uint32_t) std::atol(tmp.c_str());
  /// exon starts
  getline(line, tmp, '\t');
  exonStarts = splitStringToInt(tmp, ",");
  /// exon ends
  getline(line, tmp, '\t');
  exonEnds = splitStringToInt(tmp, ",");
  /// score
  getline(line, tmp, '\t');
  score = (int) std::atol(tmp.c_str());
  /// gene name
  getline(line, tmp, '\t');
  geneName = tmp;
  /// start status
  getline(line, tmp, '\t');
  cdsStartStat = tmp;
  /// end status
  getline(line, tmp, '\t');
  cdsEndStat = tmp;
  /// exon frames
  getline(line, tmp, '\t');
  exonFrames = tmp;
  
  /// start position is 0-based inclusive
  /// end position is 0-based exclusive, ie, [start,end)
  transcriptLength = txEnd - txStart;
  removeUTR();
}

/**
 * de-constructor
 */
RefSeqTranscript::~RefSeqTranscript()
{
  
}

/**
 * @brief remove the UTR in exon coordinates
 */
void RefSeqTranscript::removeUTR()
{
  cDNALength = 0;
  codingExonCout = 0;
  codingExonStarts.clear();
  codingExonEnds.clear();
  /// only for NM transcript
  if (cdsStart != cdsEnd)
  {
    for (int index = 0; index < exonCount; ++index)
    {
      if ((exonStarts[index] < cdsEnd) && (exonEnds[index] > cdsStart))
      {
        if ((exonStarts[index] < cdsStart) && (exonEnds[index] > cdsStart) && exonEnds[index] <= cdsEnd)
        {
          codingExonStarts.push_back(cdsStart);
          codingExonEnds.push_back(exonEnds[index]);
        }
        else if (exonStarts[index] < cdsEnd && exonEnds[index] > cdsEnd && exonStarts[index] >= cdsStart)
        {
          codingExonStarts.push_back(exonStarts[index]);
          codingExonEnds.push_back(cdsEnd);
          
        }
        else if (exonEnds[index] > cdsEnd && exonStarts[index] < cdsStart)
        {
          codingExonStarts.push_back(cdsStart);
          codingExonEnds.push_back(cdsEnd);
        }
        else
        {
          codingExonStarts.push_back(exonStarts[index]);
          codingExonEnds.push_back(exonEnds[index]);
        }
        
      }
      
    }
    
    codingExonCout = (int) codingExonStarts.size();
    
    for (int i = 0; i < codingExonStarts.size(); ++i)
    {
      cDNALength += codingExonEnds[i] - codingExonStarts[i];
      
    }
  }
}

/**
 * @brief split a string to int vector
 * @param s
 * @param delim
 * @return
 */
std::vector<uint32_t> RefSeqTranscript::splitStringToInt(const std::string &s, const std::string &delim)
{
  bool keep_empty = false;
  std::vector<uint32_t> result;
  result.empty();
  if (delim.empty())
  {
    result.push_back((const unsigned int &) std::atol(s.c_str()));
    return result;
  }
  
  std::string::const_iterator substart = s.begin(), subend;
  while (true)
  {
    subend = search(substart, s.end(), delim.begin(), delim.end());
    std::string temp(substart, subend);
    if (keep_empty || !temp.empty())
    {
      result.push_back((const unsigned int &) std::atol(temp.c_str()));
    }
    if (subend == s.end())
    {
      break;
    }
    substart = subend + delim.size();
  }
  return result;
}

/**
 * @brief Print out the exon coodinates for debuging
 * @param os
 */
void RefSeqTranscript::printOut(std::ostream &os)
{
  os << bin << "\t";
  os << transcriptID << "\t";
  os << geneName << "\t";
  os << chrom << "\t";
  os << exonCount << "\t";
  os << codingExonCout << "\t";
  os << cdsStart << "\t" << cdsEnd << "\t";
  for (auto const &pos:codingExonStarts)
    os << pos << ",";
  os << "\t";
  for (auto const &pos:codingExonEnds)
    os << pos << ",";
  os << "\t";
  os << std::endl;
  
}



/**
 * @brief read the original UCSC refGene table
 * @param filename
 * @return
 */
void readRefSeqTranscript(std::string filename, std::vector<RefSeqTranscript> &refSeqTranscripts)
{
  refSeqTranscripts.clear();
  std::ifstream in;
  /// read the original refSeq gene annotation table
  in.open(filename);
  if (!in.is_open())
  {
    std::cerr << "Error: cannot open \t" << filename << std::endl;
    exit(1);
  }
  std::string line;
  std::stringstream line2;
  std::string tmp2;
  regex e1("NR_");
  regex e2("\tJAK3\t");
  smatch m;
  while (getline(in, line, '\n'))
  {
    line2.str("");
    line2.clear();
    line2<<line;
//    if (regex_search(line,m,e2))
//    {
//      cerr<<"read JAK3!\n";
//      cerr<<line << std::endl;
//      int i = 0;
//      i++;
//      cerr<<i<<std::endl;
//      regex_search(line,m,e2);
//    }
    getline(line2, tmp2, '\t');
    getline(line2, tmp2, '\t');
    if (regex_search(tmp2,m,e1))
    {
      continue;
    }
    
    RefSeqTranscript tx(line);
//    if (tx.geneName=="JAK3")
//      cerr<<"read JAK3!\n";
    /// only consider NM transcripts
    //if (tx.transcriptID.substr(0,2) == "NM")
//    if (tx.cDNALength > 0)
    refSeqTranscripts.push_back(tx);
  }
  in.close();
  
  
}

void
find_the_longest_txpts_per_gene(vector<RefSeqTranscript> &txpts, vector<RefSeqTranscript> &longest_txpts)
{
  std::map<string, RefSeqTranscript> txpt_gene;
  std::map<string, RefSeqTranscript>::iterator txpt_gene_it;
  string gene_name = "";
  for (int j = 0; j < txpts.size(); ++j)
  {
    gene_name = txpts[j].geneName;
    txpt_gene_it = txpt_gene.find(gene_name);
    if (txpt_gene_it == txpt_gene.end())
    {
      txpt_gene[gene_name] = txpts[j];
    }
    else
    {
      if (txpts[j].cDNALength > txpt_gene_it->second.cDNALength)
      {
        txpt_gene[gene_name] = txpts[j];
      }
    }
  }
  // produce longest txpts vector
  std::map<string, RefSeqTranscript>::iterator txpt_gene_it2;
  for (txpt_gene_it2 = txpt_gene.begin(); txpt_gene_it2 != txpt_gene.end(); ++txpt_gene_it2)
  {
    longest_txpts.push_back(txpt_gene_it2->second);
  }
//  std::cout << longest_txpts.size() << std::endl;
  
//  ///temp pro to produce gene-longesttxpts
//  string temp_file= "/home/jinlf/data/gene_longest_txpts_map.txt";
//  ofstream out;
//  out.open(temp_file.c_str());
//  for (int i = 0; i < longest_txpts.size(); ++i)
//  {
//    out << longest_txpts[i].geneName << "\t" << longest_txpts[i].transcriptID<<std::endl;
//  }
//  out.close();
}

void add_cds_parts(vector<RefSeqTranscript> &txpts)
{
  for (int k = 0; k < txpts.size(); ++k)
  {
    for (int i = 0; i < txpts[k].codingExonStarts.size(); ++i)
    {
      txpts[k].codingExonParts.push_back(txpts[k].codingExonStarts[i]);
      txpts[k].codingExonParts.push_back(txpts[k].codingExonEnds[i]);
    }
  }
  
}

void find_the_longest_cds_txpt(vector<RefSeqTranscript> &txpts, RefSeqTranscript &longestCdsTxpt)
{
  int max_cds_length = 0;
  for (int i = 0; i < txpts.size(); ++i)
  {
    if (txpts[i].cDNALength >max_cds_length){
      longestCdsTxpt = txpts[i];
    }
  }
}


void find_genes_from_txpts(vector<RefSeqTranscript> &txpts_vec, set<string> &genes_vec)
{
  genes_vec.clear();
  for (int i = 0; i < txpts_vec.size(); ++i)
  {
    genes_vec.insert(txpts_vec[i].geneName);
  }
}

