/* 
 * CrossbarSHD
 * Software to simulate the CrossbarSHD algorithm
 *
 * output: Read_Name    0/16    Chromosome_name   Index_in_chromosome   Error_count
 * 
 * crossbar_shd.cpp
 * Author: Keeton Hodgson
 * Last modified: 1/27/2017
 */

#include <stdio.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>

#include "compare.h"
#include "reference.h"
#include "commandline.h"

#include "ctpl.h"

using namespace std;

void read_compare_func16(int id, string header, string read_line, string chrom_name, vector<unsigned short> * ref, int shift, int threshold);
void read_compare_func4(int id, string header, string read_line, string chrom_name, vector<unsigned char> * ref, int shift, int threshold);

void CompareRead16(char* read_file, char* reference_file, int shift, int threshold);
void CompareRead4(char* read_file, char* reference_file, int shift, int threshold);

char* sequence = 0;
int threads = 2;

int main(int argc, char** argv)
{
  if(argc < 2) {
    cerr << "Error: no arguments specified.\n";
    return 1;
  }

  char* reference_file = 0;
  char* read_file = 0;

  // Set some default parameters
  int shift = 0;
  int encoding = 1;
  int threshold = 0;
  
  if(CommandLineOptions(argc, argv, &shift, &encoding, &threshold, &reference_file, &read_file, &sequence, &threads))
    exit(0);

  //cout<<"shift: "<< shift<<endl;
  //cout<<"threshold: "<< threshold<<endl;

  if(!reference_file || !read_file) {
    cerr<< "Need to specify a reference and read file.\n";
    return 1;
  } 

  if(encoding){
    //cout<<"16 Bit encodings.\n";
    CompareRead16(read_file, reference_file, shift, threshold);
  }
  else{
    //cout<<"4 Bit encodings.\n";
    CompareRead4(read_file, reference_file, shift, threshold);
  }

  return 0;
}

/*
 * Compares each read in the read_file to the reference genome according to the shift and threshold
 * 16 bit version
 */
void CompareRead16(char* read_file, char* reference_file, int shift, int threshold) {
  vector<vector<unsigned short> * > ref(0);
  vector<string> ref_names(0);
  ref.reserve(64);
  ref_names.reserve(64);
  PrepareReference16(reference_file, &ref, &ref_names);
  ctpl::thread_pool p(threads);
  string str;
  if(sequence != 0)
    str = string(sequence);

  // Set up and read from the file with the read sequences
  string read_line;
  ifstream file(read_file);
  if(file.is_open()) {

    // Step through the read sequences
    while(getline(file, read_line)) {
      if(read_line.size() > 2 && read_line[0] == '@' && read_line[1] == 'E' && read_line[2] == 'R') {
        char c = read_line[1];
        string name = "";
        int cindex = 1;
        while(c != ' ') {
          name+=c;
          cindex++;
          c = read_line[cindex];
        }
        // s+='\t';
        if(sequence != 0){
          if(str == name)
            sequence = 0;
          else
            continue;
        }


        getline(file, read_line);
        for(unsigned int i = 0; i < ref.size(); i++){
          p.push(read_compare_func16, name, read_line, ref_names.at(i), ref.at(i), shift, threshold);
        }
      }
    }
  }
  else {
    cerr << "Unable to open reference file \n";
    exit(1);
  }
} 

/*
 * Helper function for CompareRead16. This function is passed into thread arguments
 */
void read_compare_func16(int id, string header, string read_line, string chrom_name, vector<unsigned short> * ref, int shift, int threshold) {
  int j = 0;
  // Save read sequence
  //header += read_line + '\n';

  std::vector<unsigned short> readv(0);
  std::vector<unsigned short> read_inverse(0);
  int read_size = read_line.size();
  for(; j < read_size-1; j++) {
    readv.push_back(ConvertCharacters16(read_line[j], read_line[j+1]));
    read_inverse.push_back(ConvertInverseCharacters16(read_line[read_size-(j+1)], read_line[read_size-(j+2)]));
  }

  int num_matches = 0;
  string res = compare16(ref, readv, read_inverse, shift, threshold, &num_matches, header, chrom_name);

  if(res.size() > 1){
    //if(report_total_matches)
      //cout<<res<<"TOTAL MATCHES: " << num_matches << '\n' << endl;
    //else
      cout<<res.substr(0, res.size()-1)<<endl;
  }
}

/*
 * Compares each read in the read_file to the reference genome according o the shift and threshold
 * 4 bit version
 */
void CompareRead4(char* read_file, char* reference_file, int shift, int threshold) {
  vector<vector<unsigned char>*> ref(0);
  vector<string> ref_names(0);
  ref.reserve(64);
  ref_names.reserve(64);
  PrepareReference4(reference_file, &ref, &ref_names);
  ctpl::thread_pool p(8);

  // Set up and read from the file with the read sequences
  string read_line;
  ifstream file(read_file);
  if(file.is_open()) {

    // Step through the read sequences
    while(getline(file, read_line)) {
      if(read_line.size() > 2 && read_line[0] == '@' && read_line[1] == 'E' && read_line[2] == 'R') {
        char c = read_line[1];
        string s = "";
        int cindex = 1;
        while(c != ' ') {
          s+=c;
          cindex++;
          c = read_line[cindex];
        }
        //s+='\t';

        getline(file, read_line);
        for(unsigned int i = 0; i < ref.size(); i++) {
          p.push(read_compare_func4, s, read_line, ref_names.at(i), ref.at(i), shift, threshold);
        }
      }
    }
  }
  else {
    cerr << "Unable to open reference file \n";
    exit(1);
  }
} 

/*
 * Helper function for CompareRead4. THis function is passed into threads to be run
 */
void read_compare_func4(int id, string header, string read_line, string chrom_name, vector<unsigned char> * ref, int shift, int threshold) {
  int j = 0;
  // Save read sequence
  //header += read_line + '\n';

  std::vector<unsigned char> readv(0);
  std::vector<unsigned char> read_inverse(0);
  int read_size = read_line.size();
  for(; j < read_size-1; j++) {
    readv.push_back(ConvertCharacter4(read_line[j]));
    read_inverse.push_back(ConvertCharacterInverse4(read_line[read_size-(j+1)]));
  }

  int num_matches = 0;
  string res = compare4(ref, readv, read_inverse, shift, threshold, &num_matches, header, chrom_name);

  if(res.size() > 1) {
  //if(report_total_matches)
    //cout<<res<<"TOTAL MATCHES: " << num_matches << '\n' << endl;
  //else
    cout<<res.substr(0, res.size() - 1)<<endl;
  }
}
