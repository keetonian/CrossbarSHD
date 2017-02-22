#include "compare.h"
#include <sstream>
#include <string>

using namespace std;



/**
 * Function to do the 4-bit comparison.
 * ref: pointer to reference genome (large compared to read). Format: 16-bit 1-hot encodings.
 * read: pointer to read sequence. Again, 4-bit 1-hot encodings.
 * shift: multiplex distance
 * threshold: minimum value allowed.
 */
string compare4(vector<unsigned char>* ref, vector<unsigned char> read_in, vector<unsigned char> inverse, int shift, int threshold, int* num_matches, string name, string chrom_name) {

  // Set up space for results.
  stringstream s;
  int matches = 0;

  // Set up variables
  unsigned int i = 0;
  unsigned int k = 0;
  unsigned int read_size = read_in.size();
  unsigned int ref_size = ref->size();
  vector<unsigned char> read(0);
  vector<unsigned char> inv(0);
  inv.reserve(read_size);
  read.reserve(read_size);

  // Prepare the forward and backwards reads
  for(i = 0; (unsigned)i < read_size; i++) {
    unsigned char a = read_in.at(i);
    unsigned char b = inverse.at(i);
    unsigned int j; 
    for(j= 0; j <= (unsigned)shift; j++) {
      if(i >= j){
        a = a | read_in.at(i-j);
        b = b | inverse.at(i-j);
      }
      if((i+j)<read_size){
        a = a | read_in.at(i+j);
        b = b | inverse.at(i+j);
      }
    }
    read.push_back(a);
    inv.push_back(b);
  }

  // Compare read strand with entire reference genome.
  for(k = 0; k <= ref_size - read_size; k++) {
    int error = 0;
    int error_inv = 0;

    for(i = 0; i < read_size; i++) {
      // Prepare the read parameters
      unsigned char a = read.at(i);
      unsigned char b = inv.at(i);
      unsigned char r = ref->at(k+i);

      // Do the comparison
      error += (!(a&r));
      error_inv += (!(b&r));
      if(error > threshold && error_inv > threshold) {
        break;
      }
    }

    //If the result is over or equal to the threshold
    if(error <= threshold){
      s << name << "\t0\t" << chrom_name << '\t' << k+1 << '\t' << error << '\t';
      s << "\n";
      matches += 1;
    }

    //If the result is over or equal to the threshold
    if(error_inv <= threshold){
      s << name << "\t16\t" << chrom_name << '\t' << k+1 << '\t' << error_inv << '\t';
      s << "\n";
      matches += 1;
    }
  }

  (*num_matches) = matches;

  return s.str();
}
