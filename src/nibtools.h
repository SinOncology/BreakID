#pragma once

#include<iostream>
#include<fstream>

#define MSK 0x6be93d3a

//error messages
const std::string errormsg[] = {"", "wrong format", "cannot open file", "file is not open",
                                "position beyond sequence boundary"};

class nib
{
private:
  long nBases;  //sequence length
  std::ifstream in;
  char raw[4]; //binary input
  unsigned long tmp;
  int base_low;
  int base_high;
  bool is_high;
  
  int bin2ascii(char *out, int inp) //converts binary to ascii
  {
    switch (inp & 0xff)
    {
    case 0:
      *out = 'T';
      break;
    case 1:
      *out = 'C';
      break;
    case 2:
      *out = 'A';
      break;
    case 3:
      *out = 'G';
      break;
    case 4:
      *out = 'N';
      break;
    case 8:
      *out = 'T';
      break;
    case 9:
      *out = 'C';
      break;
    case 10:
      *out = 'A';
      break;
    case 11:
      *out = 'G';
      break;
    default:
      *out = 'N';
      return 1;
    }
    return 0;
  };
  
  unsigned char ascii2bin(char inp)
  {
    switch (inp)
    {
    case 'T':
      return 0;
    
    case 'C':
      return 1;
    
    case 'A':
      return 2;
    
    case 'G':
      return 3;
    
    case 'N':
      return 4;
    
    case 't':
      return 8;
    
    case 'c':
      return 9;
    
    case 'a':
      return 10;
    
    case 'g':
      return 11;
    
    default:
      return 4;
    }
  };
public:
  int open(std::string filename);
  
  void close()
  {
    if (!in.is_open())
      in.close();
  }
  
  int getBase(char *base, unsigned long pos);
  
  int nextBase(char *base);
  
  int write(std::string seq, std::string filename);
  
  unsigned long size()
  {
    if (!in.is_open())
      return 0;
    else
      return nBases;
  }
};
