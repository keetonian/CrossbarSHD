#ifndef COMMANDLINE_H
#define COMMANDLINE_H

int CommandLineOptions(int argc, char** argv, int * shift, int * encoding, int * threshold, char ** reference_file, char ** read_file, char ** sequence, int * threads);

/*
 * Parses a file path from the command line arguments
 */
char* ParseFile(char** argv, int i);


#endif //COMMANDLINE_H
