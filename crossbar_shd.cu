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
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "compare.h"
#include "reference.h"
#include "commandline.h"

using namespace std;

void read_compare_func16(string header, string read_line, string chrom_name, vector<unsigned short> * ref, int shift, int threshold, int seed);
void read_compare_func4(int id, string header, string read_line, string chrom_name, vector<unsigned char> * ref, int shift, int threshold);

void CompareRead16(char* read_file, char* reference_file, int shift, int threshold);
void CompareRead4(char* read_file, char* reference_file, int shift, int threshold);

/**
 * Function to do the 16-bit comparison.
 * ref: pointer to reference genome (large compared to read). Format: 16-bit 1-hot encodings.
 * read: pointer to read sequence. Again, 16-bit 1-hot encodings.
 */
__global__ void compare16(const unsigned short * ref, const unsigned short * read, char * output, int numElements, int rIndex) {

  __shared__ short rd[1024];

  int i = threadIdx.x + threadIdx.y + blockIdx.x*blockDim.x;
  int ii = threadIdx.x + threadIdx.y*blockDim.x;

  rd[ii] = 0;
  if(i < numElements){
    if(read[threadIdx.x+rIndex]&ref[i])
      rd[ii] = 1;
  }

  __syncthreads();

  if(threadIdx.x % 2 == 0 && threadIdx.x + 1 < blockDim.x){
    rd[ii] += rd[ii+1];
  }
  __syncthreads();
  if(threadIdx.x % 4 == 0 && threadIdx.x + 2 < blockDim.x){
    rd[ii] += rd[ii+2];
  }
  __syncthreads();
  if(threadIdx.x % 8 == 0 && threadIdx.x + 4 < blockDim.x){
    rd[ii] += rd[ii+4];
  }
  __syncthreads();
  if(threadIdx.x % 16 == 0 && threadIdx.x + 8 < blockDim.x){
    rd[ii] += rd[ii+8];
  }
  __syncthreads();
  if(threadIdx.x % 32 == 0 && threadIdx.x + 16 < blockDim.x){
    rd[ii] += rd[ii+16];
  }
  __syncthreads();

  if(threadIdx.x == 0)
    output[i] = rd[ii];
}

char* sequence = 0;
int seed_size = 32;

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

  if(CommandLineOptions(argc, argv, &shift, &encoding, &threshold, &reference_file, &read_file, &sequence, &seed_size))
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

  // Convert the reference to 16-bit vectors
  vector<vector<unsigned short> * > ref(0);
  vector<string> ref_names(0);
  ref.reserve(64);
  ref_names.reserve(64);
  PrepareReference16(reference_file, &ref, &ref_names);


  // Get the starting read if it exists.
  string str;
  if(sequence != 0)
    str = string(sequence);

  // Do the comparisons one chromosome at a time to be more memory efficient
  for(unsigned int chromosome = 0; chromosome < ref.size(); chromosome++){

    string chrom_name = ref_names.at(chromosome);

    // Prepare sizes of elements
    int numElements=ref.at(chromosome)->size();
    size_t size = numElements*sizeof(short);
    size_t size_out = numElements*sizeof(char);
    char* hOutput = (char*)malloc(size_out);


    // Reserve space for the Reference, Output arrays in the GPU
    cudaError_t err = cudaSuccess;
    unsigned short *cRef = NULL;
    err = cudaMalloc((void **) &cRef, size);
    if(err != cudaSuccess){
      cerr << "Failed to allocate device vector ref (error code " << cudaGetErrorString(err) << ")!" << endl;
      exit(1);
    }

    char* cOutput = NULL;
    err = cudaMalloc((void **) &cOutput, size_out);
    if(err != cudaSuccess) {
      cerr << "CudaFailure on cOutput" << cudaGetErrorString(err) << endl;
      exit(1);
    }

    // Copy the reference chromosome into the GPU
    err = cudaMemcpy(cRef, (unsigned short*)(&(ref.at(chromosome)->at(0))), size, cudaMemcpyHostToDevice);
    if(err != cudaSuccess) {
      cerr << "Cuda Failure on copying cRef" << endl;
      exit(1);
    }

    // Prepare the thread dimensions
    dim3 threadsPerBlock(seed_size, seed_size);
    int blocksPerGrid= (numElements + seed_size)/(seed_size);

    // Extract the reads (name and sequence) from the read file
    string read_line;
    ifstream file(read_file);
    if(file.is_open()) {

      // Step through the read sequences
      while(getline(file, read_line)) {
	
	// Extract the name
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

	  // Get the sequence
	  getline(file, read_line);

	  // Prepare the sequences to be used by the GPU
	  int read_size = read_line.size()-1;
	  vector<unsigned short> readv_temp, read_inverse_temp;

	  // Allocate space for the arrays
	  unsigned short * readv, * read_inverse;
	  readv = (unsigned short *)malloc((read_size)*sizeof(short));
	  read_inverse = (unsigned short *)malloc((read_size)*sizeof(short));
	  size_t size_read = (read_size)*sizeof(short);
	  for(int j = 0; j < read_size; j++) {
	    readv_temp.push_back(ConvertCharacters16(read_line[j], read_line[j+1]));
	    read_inverse_temp.push_back(ConvertInverseCharacters16(read_line[read_size-(j)], read_line[read_size-(j+1)]));
	  }

	  // Copy the data into the arrays, using appropriate shift distance
	  for(int i = 0; i < readv_temp.size(); i++){
	    unsigned short a = readv_temp.at(i);
	    unsigned short b = read_inverse_temp.at(i);
	    for(int j = 0; j <= shift; j++){
	      if(i>=j){
		a = a | readv_temp.at(i-j);
		b = b | read_inverse_temp.at(i-j);
	      }
	      if((i+j)<readv_temp.size()){
		a = a | readv_temp.at(i+j);
		b = b | read_inverse_temp.at(i+j);
	      }
	    }
	    readv[i] = a;
	    read_inverse[i] = b;
	  }

	  // Allocate space on GPU for the forward/reverse inverse reads
	  unsigned short *cRead = NULL;
	  err = cudaMalloc((void **) &cRead, size_read);
	  if(err != cudaSuccess) {
	    cerr << "Cuda failure on cRead" << endl;
	    exit(1);
	  }

	  unsigned short *cInverse = NULL;
	  err = cudaMalloc((void **) &cInverse, size_read);
	  if(err != cudaSuccess) {
	    cerr << "Cuda failure on cInverse" << endl;
	    exit(1);
	  }

	  err = cudaMemcpy(cRead, readv, size_read, cudaMemcpyHostToDevice);
	  if(err != cudaSuccess) {
	    cerr << "Cuda failure on copying cRead " << cudaGetErrorString(err) << endl;
	    exit(1);
	  }

	  err = cudaMemcpy(cInverse, read_inverse, size_read, cudaMemcpyHostToDevice);
	  if(err != cudaSuccess) {
	    cerr << "Cuda failure on copying cInverse " << cudaGetErrorString(err) << endl;
	    exit(1);
	  }

	  // Do the calculations for each seed
	  for(int i = 0; i < (read_size)/seed_size; i++){

	    compare16<<<blocksPerGrid,threadsPerBlock>>>(cRef, cRead, cOutput, numElements, i*seed_size);
	    err = cudaGetLastError();
	    if(err != cudaSuccess){
	      cout << "Failure on cuda computation: " << cudaGetErrorString(err) << endl;
	      exit(1);
	    }

	    err = cudaMemcpy(hOutput, cOutput, size_out, cudaMemcpyDeviceToHost);
	    if(err != cudaSuccess){
	      cerr << "Failure on copying output" << cudaGetErrorString(err) << endl;
	      exit(1);
	    }

	    // Print the results
	    for(int ii = 0; ii < numElements; ii++)
	      if(hOutput[ii] >= seed_size-threshold)
		cout << name << '\t' << seed_size << '\t' << i << '\t' << "0" << '\t' << chrom_name << '\t' << ii << '\t' << (int)hOutput[ii] << endl;

	  }

	  // Do the calculations for each rinverse seed
	  for(int i = 0; i < (read_size)/seed_size; i++){

	    compare16<<<blocksPerGrid,threadsPerBlock>>>(cRef, cInverse, cOutput, numElements, i*seed_size);
	    err = cudaGetLastError();
	    if(err != cudaSuccess){
	      cout << "Failure on cuda computation: " << cudaGetErrorString(err) << endl;
	      exit(1);
	    }

	    err = cudaMemcpy(hOutput, cOutput, size_out, cudaMemcpyDeviceToHost);
	    if(err != cudaSuccess){
	      cerr << "Failure on copying output" << cudaGetErrorString(err) << endl;
	      exit(1);
	    }

	    // Print the results
	    for(int ii = 0; ii < numElements; ii++)
	      if(hOutput[ii] >= seed_size-threshold)
		cout << name << '\t' << seed_size << '\t' << i << '\t' << "16" << '\t' << chrom_name << '\t' << ii << '\t' << (int)hOutput[ii] << endl;
	  }

	  err = cudaFree(cRead);
	  if(err != cudaSuccess) {
	    cerr << "Failure freeing read" << endl;
	    exit(1);
	  }
	  err = cudaFree(cInverse);
	  if(err != cudaSuccess) {
	    cerr << "failure freeing cinverse" << endl;
	    exit(1);
	  }
	  free(readv);
	  free(read_inverse);

	}
      }
    }
    else {
      cerr << "Unable to open reference file \n";
      exit(1);
    }

    cout << "HERE" << endl;

    // Free all GPU memory
    err = cudaFree(cRef);
    if(err != cudaSuccess) {
      cerr << "Failure freeing reference" << endl;
      exit(1);
    }

    err = cudaFree(cOutput);
    if(err != cudaSuccess) {
      cerr << "Failure freeing output" << endl;
      exit(1);
    }

    free(hOutput);
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
  //ctpl::thread_pool p(8);

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
	  //      p.push(read_compare_func4, s, read_line, ref_names.at(i), ref.at(i), shift, threshold);
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
