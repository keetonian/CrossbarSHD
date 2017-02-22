#include "commandline.h"
#include <iostream>

using namespace std;

int CommandLineOptions(int argc, char** argv, int * shift, int * encoding, int * threshold, char ** reference_file, char ** read_file, char ** sequence, int * threads){
  for(int i = 0; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
        case 'g': //reference genome
          {
            *reference_file = ParseFile(argv, i);
            if(*reference_file == 0) {
              cerr << "Error parsing reference file.\n";
              return 1;
            }
            // Debugging output, remove later.
            // cout << "Reference: " << reference_file << "\n";
          }
          break;
        case '4': (*encoding) = 0;
                  break;
        case '1': (*encoding) = 1;
                  break;
        //case 'm': report_total_matches = 1;
                  //break;
        case 'h': //help
                  {
                    cout << "This program compares a read sequence with a reference genome.\n";
                    cout << "\t-g: reference genome path\n";
                    cout << "\t-r: read sequence path\n";
                    cout << "\t-h: display this information\n";
                    cout << "\t-t: tolerated error threshold\n";
                    cout << "\t\t0: only exact matches (no errors tolerated)\n";
                    cout << "\t-s: shift distance\n";
                    cout << "\t-4: use 4-bit encodings\n";
                    cout << "\t-16: use 16-bit encodings\n";
                    cout << "\t-n: number of threads to use (default = 2)\n";
                    cout << "\t-p: name of the read sequence to start from \n\t\t(default: none, start from beginning of file)\n";
                    //cout << "\t-m: turn on total match reporting\n";
                    //cout << "\t-o1: output grouped by sequence name\n";
                    //cout << "\t-o2: output grouped by chromosome name\n";
                    // Using the help flag will only print this help dialog.
                    return 0;
                  }
                  break;
        case 'r': //read file
                  {
                    *read_file = ParseFile(argv, i);
                    if(*read_file == 0) {
                      cerr << "Error parsing read file.\n";
                      return 1;
                    }
                  }
                  break;
        case 't': //threshold (read size - threshold) 0 means no errors tolerated
                  {
                    if(argv[i][2] != 0)
                      (*threshold) = atoi(argv[i]+2);
                    else
                      (*threshold) = atoi(argv[i+1]);
                  }
                  break;
        case 's': //shift distance
                  {
                    if(argv[i][2] != 0)
                      (*shift) = atoi(argv[i]+2);
                    else
                      (*shift) = atoi(argv[i+1]);
                  }
                  break;
        case 'p': // starting sequence
                  {
                    *sequence = ParseFile(argv, i);
                    if(*sequence == 0) {
                      cerr << "Error: specify a sequence after -p" << endl;
                      return 1;
                    }
                  }
                  break;
        case 'n': // number threads
                  {
                    if(argv[i][2] != 0)
                      *threads = atoi(argv[i]+2);
                    else
                      *threads = atoi(argv[i+1]);
                  }
                  break;
        default: //unknown flag
                  {
                    cerr << "Error: illegal flag: " << argv[i] << "\n";
                    cerr << "Use -h to print help dialog.\n";
                    return 1;
                  }
                  break;
      }
    }
  }

  return 0;
  //cout<<"shift: "<< shift<<endl;
  //cout<<"threshold: "<< threshold<<endl;
}

/*
 * Parses a file path from the command line arguments
 */
char* ParseFile(char** argv, int i) {
  char* file_name = 0;
  if(argv[i][2] == 0) {
    //Means that the file should be in the next argument.
    file_name = argv[i+1];
    if(file_name[0] == '-') {
      return 0;
    }
  }
  else {
    //Means that the file should be part of the current argument.
    file_name = argv[i] + 2;
  }
  return file_name;
}

