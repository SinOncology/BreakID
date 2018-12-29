/**
 * @file CigarRoller.cc
 * @brief Cigar Roller class
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "CigarRoller.h"


CigarRoller &CigarRoller::operator+=(CigarRoller &rhs)
{
  std::vector<CigarOperator>::iterator i;
  for (i = rhs.cigarOperations.begin(); i != rhs.cigarOperations.end(); i++)
  {
    (*this) += *i;
  }
  return *this;
}


//
// Append a new operator at the end of the sequence.
//
CigarRoller &CigarRoller::operator+=(const CigarOperator &rhs)
{
  // Adding to the cigar, so the query & reference indexes would be
  // incomplete, so just clear them.
  clearQueryAndReferenceIndexes();
  
  if (rhs.count == 0)
  {
    // nothing to do
  }
  else if (cigarOperations.empty() || cigarOperations.back() != rhs)
  {
    cigarOperations.push_back(rhs);
  }
  else
  {
    // last stored operation is the same as the new one, so just add it in
    cigarOperations.back().count += rhs.count;
  }
  return *this;
}


CigarRoller &CigarRoller::operator=(CigarRoller &rhs)
{
  clear();
  
  (*this) += rhs;
  
  return *this;
}


//
void CigarRoller::Add(Operation operation, int count)
{
  CigarOperator rhs(operation, (uint32_t) count);
  (*this) += rhs;
}


void CigarRoller::Add(char operation, int count)
{
  switch (operation)
  {
  case 0:
  case 'M':
    Add(match, count);
    break;
  case 1:
  case 'I':
    Add(insert, count);
    break;
  case 2:
  case 'D':
    Add(del, count);
    break;
  case 3:
  case 'N':
    Add(skip, count);
    break;
  case 4:
  case 'S':
    Add(softClip, count);
    break;
  case 5:
  case 'H':
    Add(hardClip, count);
    break;
  case 6:
  case 'P':
    Add(pad, count);
    break;
  case 7:
  case '=':
    Add(match, count);
    break;
  case 8:
  case 'X':
    Add(match, count);
    break;
  default:
    // Hmmm... what to do?
    std::cerr << "ERROR "
              << "(" << __FILE__ << ":" << __LINE__ << "): "
              << "Parsing CIGAR - invalid character found "
              << "with parameter " << operation << " and " << count
              << std::endl;
    break;
  }
}


void CigarRoller::Add(const char *cigarString)
{
  int operationCount = 0;
  while (*cigarString)
  {
    if (isdigit(*cigarString))
    {
      char *endPtr;
      operationCount = (int) strtol((char *) cigarString, &endPtr, 10);
      cigarString = endPtr;
    }
    else
    {
      Add(*cigarString, operationCount);
      cigarString++;
    }
  }
}


bool CigarRoller::Remove(int index)
{
  if ((index < 0) || ((unsigned int) index >= cigarOperations.size()))
  {
    // can't remove, out of range, return false.
    return (false);
  }
  cigarOperations.erase(cigarOperations.begin() + index);
  // Modifying the cigar, so the query & reference indexes are out of date,
  // so clear them.
  clearQueryAndReferenceIndexes();
  return (true);
}


bool CigarRoller::IncrementCount(int index, int increment)
{
  if ((index < 0) || ((unsigned int) index >= cigarOperations.size()))
  {
    // can't update, out of range, return false.
    return (false);
  }
  cigarOperations[index].count += increment;
  
  // Modifying the cigar, so the query & reference indexes are out of date,
  // so clear them.
  clearQueryAndReferenceIndexes();
  return (true);
}


bool CigarRoller::Update(int index, Operation op, int count)
{
  if ((index < 0) || ((unsigned int) index >= cigarOperations.size()))
  {
    // can't update, out of range, return false.
    return (false);
  }
  cigarOperations[index].operation = op;
  cigarOperations[index].count = (uint32_t) count;
  
  // Modifying the cigar, so the query & reference indexes are out of date,
  // so clear them.
  clearQueryAndReferenceIndexes();
  return (true);
}


void CigarRoller::Set(const char *cigarString)
{
  clear();
  Add(cigarString);
}


void CigarRoller::Set(const uint32_t *cigarBuffer, uint16_t bufferLen)
{
  clear();
  
  // Parse the buffer.
  for (int i = 0; i < bufferLen; i++)
  {
    int opLen = cigarBuffer[i] >> 4;
    
    Add(cigarBuffer[i] & 0xF, opLen);
  }
}


//
// when we examine CIGAR strings, we need to know how
// many cumulative insert and delete positions there are
// so that we can adjust the read location appropriately.
//
// Here, we iterate over the vector of CIGAR operations,
// summaring the count for each insert or delete (insert
// increases the offset, delete decreases it).
//
// The use case for this is when we have a genome match
// position based on an index word other than the first one,
// and there is also a insert or delete between the beginning
// of the read and the index word.  We can't simply report
// the match position without taking into account the indels,
// otherwise we'll be off by N where N is the sum of this
// indel count.
//
// DEPRECATED - do not use.  There are better ways to accomplish that by using
// read lengths, reference lengths, span of the read, etc.
int CigarRoller::getMatchPositionOffset()
{
  int offset = 0;
  std::vector<CigarOperator>::iterator i;
  
  for (i = cigarOperations.begin(); i != cigarOperations.end(); i++)
  {
    switch (i->operation)
    {
    case insert:
      offset += i->count;
      break;
    case del:
      offset -= i->count;
      break;
      // TODO anything for case skip:????
    default:
      break;
    }
  }
  return offset;
}


//
// Get the string reprentation of the Cigar operations in this object.
// Caller must delete the returned value.
//
const char *CigarRoller::getString()
{
  // NB: the exact size of the string is not important, it just needs to be guaranteed
  // larger than the largest number of characters we could put into it.
  
  // we do not explicitly manage memory usage, and we expect when program exits, the memory used here will be freed
  static char *ret = NULL;
  static uint32_t retSize = 0;
  
  if (ret == NULL)
  {
    retSize = (uint32_t) cigarOperations.size() * 12 + 1;  // 12 == a magic number -> > 1 + log base 10 of MAXINT
    ret = (char *) malloc(sizeof(char) * retSize);
    assert(ret != NULL);
    
  }
  else
  {
    // currently, ret pointer has enough memory to use
    if (retSize > cigarOperations.size() * 12 + 1)
    {
    }
    else
    {
      retSize = (uint32_t) cigarOperations.size() * 12 + 1;
      free(ret);
      ret = (char *) malloc(sizeof(char) * retSize);
    }
    assert(ret != NULL);
  }
  
  char *ptr = ret;
  char buf[12];   // > 1 + log base 10 of MAXINT
  
  std::vector<CigarOperator>::iterator iterator;
  
  // Progressively append the character representations of the operations to
  // the cigar string we allocated above.
  
  *ptr = '\0';    // clear result string
  for (iterator = cigarOperations.begin(); iterator != cigarOperations.end(); iterator++)
  {
    sprintf(buf, "%d%c", (*iterator).count, (*iterator).getChar());
    strcat(ptr, buf);
    while (*ptr)
    {
      ptr++;    // limit the cost of strcat above
    }
  }
  return ret;
}


void CigarRoller::clear()
{
  // Clearing the cigar, so the query & reference indexes are out of
  // date, so clear them.
  clearQueryAndReferenceIndexes();
  cigarOperations.clear();
}

uint32_t CigarRoller::getAlignmentEnd(uint32_t start_pos)
{
//  getRefPosition(,start_pos);
  
  return start_pos + getExpectedReferenceBaseCount() - 1;
}

bool CigarRoller::is_complementary_cigar(const std::string c2, int error_num)
{
  bool match = false;
  regex regex1("([0-9]+[MS]){2}");
  std::smatch sm;
  string c1;
  getCigarString(c1);
  int c1_m, c1_s, c2_m, c2_s;
  CigarRoller c2_roller;
  c2_roller.Set(c2.c_str());
  if (regex_match(c1, regex1) && regex_match(c2, regex1))
  {
    c1_m = getNumMatches();
    c2_m = c2_roller.getNumMatches();
    c1_s = getNumBeginClips() + getNumEndClips();
    c2_s = c2_roller.getNumEndClips() + c2_roller.getNumBeginClips();
    if ((c1_m <= c2_s + error_num && c1_m >= c2_s - error_num) && (c1_m + c1_s == c2_m + c2_s))
    {
      match = true;
    }
  }
  
  return match;
}

