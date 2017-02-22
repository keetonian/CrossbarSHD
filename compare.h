#ifndef COMPARE_H
#define COMPARE_H
#include <vector>
#include <string>

using namespace std;

/**
 * Function to do the 16-bit comparison.
 * ref: pointer to reference genome (large compared to read). Format: 16-bit 1-hot encodings.
 * read: pointer to read sequence. Again, 16-bit 1-hot encodings.
 * shift: multiplex distance
 * threshold: minimum value allowed.
 */
string compare16(vector<unsigned short>* ref, vector<unsigned short> read_in, vector<unsigned short> inverse, int shift, int threshold, int* num_matches, string name, string chrom_name);


/**
 * Function to do the 4-bit comparison.
 * ref: pointer to reference genome (large compared to read). Format: 16-bit 1-hot encodings.
 * read: pointer to read sequence. Again, 4-bit 1-hot encodings.
 * shift: multiplex distance
 * threshold: minimum value allowed.
 */
string compare4(vector<unsigned char>* ref, vector<unsigned char> read_in, vector<unsigned char> inverse, int shift, int threshold, int* num_matches, string name, string chrom_name);

#endif //COMPARE_H
