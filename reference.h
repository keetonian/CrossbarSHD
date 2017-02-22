#ifndef REFERENCE_H
#define REFERENCE_H

#include <vector>
#include <string>

using namespace std;

/*
 * Prepares the 16 bit reference
 * */
void PrepareReference16(char* ref_file, vector<vector<unsigned short> * > * ref, vector<string> * ref_names);

/*
 * Converts 2 characters to their 16-bit encodings
 * */
unsigned short ConvertCharacters16(char char1, char char2);

/*
 * Converts 2 characters to their 16-bit encodings
 * */
unsigned short ConvertInverseCharacters16(char char1, char char2);

/*
 * Converts the reference into 4 bit codes
 */
void PrepareReference4(char* ref_file, vector<vector<unsigned char> *> * ref, vector<string> *ref_names);

/* 
 * Converts 1 character to a 4-bit vector
 * */
unsigned char ConvertCharacter4(char char1);

/* 
 * Converts 1 character to its inverse 4-bit vector
 * */
unsigned char ConvertCharacterInverse4(char char1);

#endif
