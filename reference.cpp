#include "reference.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

/*
 * Opens the reference file, converts it to 16-bit encodings.
 * This method is very memory-hungry. It attemps to load the entire reference.
 * AA: 0x8000 
 * AT: 0x4000
 * AC: 0x2000
 * AG: 0x1000
 * TA: 0x0800
 * TT: 0x0400
 * TC: 0x0200
 * TG: 0x0100
 * CA: 0x0080
 * CT: 0x0040
 * CC: 0x0020
 * CG: 0x0010
 * GA: 0x0008
 * GT: 0x0004
 * GC: 0x0002
 * GG: 0x0001
 */
void PrepareReference16(char* ref_file, vector<vector<unsigned short> * > * ref, vector<string> * ref_names){
  string line;
  ifstream reference(ref_file);
  if(reference.is_open()) {
    char previous = 0;

    vector<unsigned short> *chromosome;// = new vector<unsigned short>();
    //cout << "parsing reference" << endl;
    while(getline(reference, line)) {
      if(line.size() ==0)
	continue;
      if (line[0] != '>'){
	unsigned int i = 0;
	if(previous)
	  chromosome->push_back(ConvertCharacters16(previous, line[0]));
	for(; i < line.size()-1; i++) {
	  chromosome->push_back(ConvertCharacters16(line[i], line[i+1]));
	}
	previous = line[line.size()-1];
      }
      else {
	string name = "";
	unsigned int j = 1;
	while(line[j] != ' ' && line[j] != '\t' && j < line.size()){
	  name+=line[j];
	  j++;
	}
	chromosome = new vector<unsigned short>();
	ref->push_back(chromosome);
	ref_names->push_back(name);
	// Store Chromosome data here
	// Note: A strand cannot match from the end of one chromosome to the beginning of the next, or can it?
      }
    }

    reference.close();
  }
  else {
    cerr << "Unable to open reference file\n";
    //exit(1);
  }
}

/*
 * Converts 2 characters to their 16-bit encodings
 * */
unsigned short ConvertCharacters16(char char1, char char2) {
  char1 = toupper(char1);
  char2 = toupper(char2);
  if(char1 == 'A') {
    switch(char2) {
      case 'A': return 0x8000;
      case 'T': return 0x4000;
      case 'C': return 0x2000;
      case 'G': return 0x1000;
      default: return 0;
    }
  }
  else if(char1 == 'T') {
    switch(char2) {
      case 'A': return 0x0800;
      case 'T': return 0x0400;
      case 'C': return 0x0200;
      case 'G': return 0x0100;
      default: return 0;
    }
  }
  else if(char1 == 'C') {
    switch(char2) {
      case 'A': return 0x0080;
      case 'T': return 0x0040;
      case 'C': return 0x0020;
      case 'G': return 0x0010;
      default: return 0;
    }
  }
  else if(char1 == 'G') {
    switch(char2) {
      case 'A': return 0x0008;
      case 'T': return 0x0004;
      case 'C': return 0x0002;
      case 'G': return 0x0001;
      default: return 0;
    }
  }
  else
    return 0;
}

/*
 * Converts 2 characters to their 16-bit encodings
 * */
unsigned short ConvertInverseCharacters16(char char1, char char2) {
  char1 = toupper(char1);
  char2 = toupper(char2);
  if(char1 == 'T') {
    switch(char2) {
      case 'T': return 0x8000;
      case 'A': return 0x4000;
      case 'G': return 0x2000;
      case 'C': return 0x1000;
      default: return 0;
    }
  }
  else if(char1 == 'A') {
    switch(char2) {
      case 'T': return 0x0800;
      case 'A': return 0x0400;
      case 'G': return 0x0200;
      case 'C': return 0x0100;
      default: return 0;
    }
  }
  else if(char1 == 'G') {
    switch(char2) {
      case 'T': return 0x0080;
      case 'A': return 0x0040;
      case 'G': return 0x0020;
      case 'C': return 0x0010;
      default: return 0;
    }
  }
  else if(char1 == 'C') {
    switch(char2) {
      case 'T': return 0x0008;
      case 'A': return 0x0004;
      case 'G': return 0x0002;
      case 'C': return 0x0001;
      default: return 0;
    }
  }
  else
    return 0;
}

/*
 * Converts the reference into 4 bit codes
 */
void PrepareReference4(char* ref_file, vector<vector<unsigned char> *> * ref, vector<string> *ref_names){
  string line;
  ifstream reference(ref_file);
  if(reference.is_open()) {

    vector<unsigned char> *chromosome;

    // What if we tried buffers instead of lines?
    while(getline(reference, line)) {
      if(line.size() ==0)
	continue;
      if (line[0] != '>'){
	unsigned int i = 0;
	for(; i < line.size(); i++) {
	  chromosome->push_back(ConvertCharacter4(line[i]));
	}
      }
      else {
	string name = "";
	int j = 1;
	while(line[j] != ' ' && line[j] != '\t') {
	  name += line[j];
	  j++;
	}
	chromosome = new vector<unsigned char>();
	ref->push_back(chromosome);
	ref_names->push_back(name);
      }
    }

    reference.close();
  }
  else {
    cerr << "Unable to open reference file\n";
    //exit(1);
  }
}


unsigned char ConvertCharacter4(char char1) {
  char1 = toupper(char1);
  switch(char1) {
    case 'A': return 0x8;
    case 'T': return 0x4;
    case 'C': return 0x2;
    case 'G': return 0x1;
    default: return 0;
  }
}

unsigned char ConvertCharacterInverse4(char char1) {
  char1 = toupper(char1);
  switch(char1) {
    case 'T': return 0x8;
    case 'A': return 0x4;
    case 'G': return 0x2;
    case 'C': return 0x1;
    default: return 0;
  }
}

