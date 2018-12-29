
#include "BamAlignment.h"



BamAlignment::BamAlignment()
  : isDuplicatedBamRecord(false)
{
  record = NULL;
  myAlignmentLength = 0;
}

/**
 * @brief constructor with  pointer to a bam record
 */

BamAlignment::BamAlignment(bam1_t *b)
  : isDuplicatedBamRecord(false)
{
  record = b;
  /// parse the cigar string and
  /// set the alignment length that this read spans in the reference
  parseCigarBinary();
}

/**
 * @brief constructor with a bam record and an indicator whether to duplicate the alignment or not
 */

BamAlignment::BamAlignment(bam1_t *b, bool duplicateBamRecord)
{
  if (duplicateBamRecord)
  {
    /// duplicate the bam record
    record = bam_dup1(b);
    isDuplicatedBamRecord = true;
  }
  else
  {
    /// just pointer to the bam record
    record = b;
    isDuplicatedBamRecord = false;
  }
  /// parse the cigar string and
  /// set the alignment length that this read spans in the reference
  parseCigarBinary();
  
}

/**
 * @brief De-constructor, clear the containers
 */
BamAlignment::~BamAlignment()
{
  
  if (isDuplicatedBamRecord)
  {
    bam_destroy1(record);
    record = NULL;
  }
  else
  {
    record = NULL;
  }
}

/**
 * @brief initiate from a bam1_t record
 * @param b
 */

void BamAlignment::readBamRecord(bam1_t *b)
{
  reset();
  record = b;
  isDuplicatedBamRecord = false;
  parseCigarBinary();
  
  
}


void BamAlignment::reset()
{
  record = NULL;
  myAlignmentLength = 0;
  isDuplicatedBamRecord = false;
}




/**
 * @brief Get 0-based inclusive leftmost position
 * @return 0-based inclusive leftmost position (start position)
 */
int32_t BamAlignment::get0BasedAlignmentStart()
{
  return record->core.pos;
}


/**
 * @brief Get 0-based inclusive rightmost position
 * @return 0-based inclusive rightmost position
 */
int32_t BamAlignment::get0BasedAlignmentEnd()
{
  /// if the alignment length is 0, return the start position
  if (myAlignmentLength == 0)
  {
    return record->core.pos;
  }
  return (record->core.pos + myAlignmentLength - 1);
  
}

/**
 *
 * @return The read name
 */
std::string BamAlignment::getReadName()
{
  return (std::string) bam_get_qname(record);
}



std::string BamAlignment::originCigar()
{
  std::string originCigar = "";
  uint8_t *bc = bam_aux_get(record, "OC"); //pointer to the barcode
  if (bc != NULL)
  {
    originCigar = ((char *) bam_aux_get(record, "OC"));
    /// ignore "Z"
    originCigar = originCigar.substr(1);
  }
  return originCigar;
}




/**
 * @brief get the cigar string
 * @return cigar string
 */
std::string BamAlignment::getCigarString()
{
  return myCigarRoller.getString();
}



/**
 * @brief Parse the cigar string and store the information in myCigarRoller
 */
void BamAlignment::parseCigarBinary()
{
  uint32_t *cigarPtr = bam_get_cigar(record);
  myCigarRoller.Set(cigarPtr, (uint16_t) record->core.n_cigar);
  //myCigarRoller.getCigarString(myCigar);
  myAlignmentLength = myCigarRoller.getExpectedReferenceBaseCount();
  myUnclippedStartOffset = myCigarRoller.getNumBeginClips();
  myUnclippedEndOffset = myCigarRoller.getNumEndClips();
  
}



long BamAlignment::getAlignmentStart()
{
  return get0BasedAlignmentStart() + 1;
}

long BamAlignment::getAlignmentEnd()
{
  return get0BasedAlignmentEnd() + 1;
}

std::string BamAlignment::getChrName()
{
  return chromID2ChrName(getReferenceID());
}


long BamAlignment::getAlignmentEndContainSoft()
{
  return get0BasedUnclippedEnd() + 1;
}





std::string BamAlignment::saTag()
{
  std::string  sa_tag;
  uint8_t *bc = bam_aux_get(record, "SA"); //pointer to the barcode
  if (bc != NULL)
  {
    sa_tag = ((char *)bc);
    /// ignore "Z"
    sa_tag = sa_tag.substr(1);
  }
  return sa_tag;
}


