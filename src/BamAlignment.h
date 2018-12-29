#pragma once
#include "sam.h"
#include <string>
#include <vector>
#include <iostream>
#include <htslib/sam.h>
#include "CigarRoller.h"
#include "Cigar.h"

#include "util_bam.h"
#include "util_bed.h"

typedef struct
{
  vector<int32_t> insertLeftRefPos;
  vector<int32_t> insertRightRefPos;
  vector<int> insertQueryIndex;
  vector<char> insertBase;
} insert_info;

/// @brief An extend bam class to store more information about an alignment
class BamAlignment {
public:

  /// constructors
  BamAlignment();
  BamAlignment(bam1_t *b);
  BamAlignment(bam1_t *b,bool duplicateBamRecord);
  /// de-constructor
  ~BamAlignment ();
  /// initiate an alignment
  void readBamRecord(bam1_t *b);
  
  /// duplicate an alignment

  void reset();

  
  /// Return the reference position associated with the specified query index
  /// or INDEX_NA based on this cigar and the specified queryStartPos which
  /// is the leftmost mapping position of the first matching base in the
  /// query.

  /// Returns the 0-based inclusive leftmost position of the alignment
  int32_t get0BasedAlignmentStart();
  
  /// Returns the 0-based inclusive rightmost position of the alignment
  int32_t get0BasedAlignmentEnd();
  
  /// get the read name
  std::string getReadName();
  
  /// get cigar string
  std::string getCigarString();
  

  
  /// get the barcode

  std::string originCigar();
 
  std::string saTag();

  
  /// Parse the cigar string to calculate the cigar length and alignment end
  void parseCigarBinary();


  int32_t  getReferenceID() {return record->core.tid;}

  int32_t  get0BasedUnclippedEnd()
  {
    // myUnclippedEndOffset will be set by get0BasedAlignmentEnd if the
    // cigar has not yet been parsed, so no need to check it here.
    return (get0BasedAlignmentEnd() + myUnclippedEndOffset);
  }


  std::string getChrName();
  

  long getAlignmentStart();
  
  long getAlignmentEnd();
  

  long getAlignmentEndContainSoft();


private:
  /// a pointer to the bam record containing information about an alignment in the bam file
  bam1_t *record;
  bool isDuplicatedBamRecord;
  /// The length of the alignment that the read spans in the reference genome.
  int myAlignmentLength;
  //std::string myCigar;
  CigarRoller myCigarRoller;
  // Unclipped alignment positions.
  int32_t myUnclippedStartOffset;
  int32_t myUnclippedEndOffset;
  insert_info myInsertInfo;
  
  
};


