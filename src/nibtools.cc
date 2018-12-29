#include <stdlib.h>
#include<iostream>
#include<fstream>

#include "nibtools.h"

int nib::open(std::string filename)
{
  unsigned long i;
  //check if file is already open
  if (in.is_open())
    in.close();
  
  in.open(filename.c_str(), std::ios::binary);
  if (!in.is_open())
    return 2;
  
  //check if the signature is 0x6be93d3a
  in.read(raw, 4);
  tmp = (raw[0] & 0xff) + ((raw[1] & 0xff) << 8) + ((raw[2] & 0xff) << 16) + ((raw[3] & 0xff) << 24);
  if (tmp != MSK)
    return 1;
  
  //get number of bases
  in.read(raw, 4);
  nBases = (raw[0] & 0xff) + ((raw[1] & 0xff) << 8) + ((raw[2] & 0xff) << 16) + ((raw[3] & 0xff) << 24);
  //check if the number of bases is the size of the file
  in.seekg(0, std::ios::end);
  i = in.tellg() * 2 - 16;
  if (!(nBases == i || nBases + 1 == i))
    return 1;
  in.seekg(8);
  is_high = 0;
  
  return 0;
}

int nib::getBase(char *base, unsigned long pos)
{
  int b1, b2, stat;
  unsigned long filePos;
  
  if (!in.is_open())
    return 3;
  if (pos >= nBases)
    return 4;
  
  //determine file position
  filePos = 8 + pos / 2;
  in.seekg(filePos);
  in.read(raw, 1);
  b1 = (raw[0] & 0xff) >> 4;
  b2 = (raw[0] & 0x0f);
  
  if (pos - (pos / 2) * 2 == 0)
    stat = bin2ascii(base, b1);
  else
    stat = bin2ascii(base, b2);
  
  if (stat != 0)
    return stat;
  
  return 0;
}

int nib::nextBase(char *base)
{
  if (is_high == 0)
  {
    in.read(raw, 1);
    base_low = (raw[0] & 0xff) >> 4;
    base_high = (raw[0] & 0x0f);
    bin2ascii(base, base_low);
    is_high = 1;
  }
  else
  {
    bin2ascii(base, base_high);
    is_high = 0;
  }
  return 0;
}
