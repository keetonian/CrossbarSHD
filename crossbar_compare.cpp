#include <unordered_map>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

void read_input_data(char* filename, unordered_map<string, unordered_map<string, vector<long>* >* > *matches);
void read_stdin(unordered_map<string, unordered_map<string, vector<long>* >* > *matches);

void compare_read_matches(unordered_map<string, unordered_map<string, vector<long>* >* >* cxdata, unordered_map<string, unordered_map<string, vector<long>* >* >* aldata);

void compare_chromosome_indexes(string read_name, unordered_map<string, vector<long>* > * cxdata, unordered_map<string, vector<long>* > * aldata);

bool search_vector(vector<long> * source, long x); 
long count_false_positives(unordered_map<string, vector<long>* > * cxdata); 
  
int main(int argc, char** argv) {
  if(argc < 2) {
    cout<<"Usage: ./cxcompare [file1] [file2]" << endl;
    cout<<"file1: CrossbarSHD output\nfile2: Alignment data (in SAM format)" << endl;
    cout<<"Alternatively file2 can be read from stdin"<<endl;
    return 0;
  }

  // Map: Map<read name, Map<chromosome name, indexes in chromosome> >
  unordered_map<string, unordered_map<string, vector<long>* >* > crossbar_data, alignment_data;
  char* file1 = argv[1];
  char* file2 = 0;
  if(argc>2)
    file2 = argv[2];

  read_input_data(file1, &crossbar_data);

  if(argc>2)
    read_input_data(file2, &alignment_data);
  else
    read_stdin(&alignment_data);

  compare_read_matches(&crossbar_data, &alignment_data);

  return 0;
}


void read_input_data(char* filename, unordered_map<string, unordered_map<string, vector<long>* >* > *matches){
  string line;
  ifstream file(filename);
  if(file.is_open()) {
    while(getline(file, line)) {
      if(line.size() == 0)
        continue;

      unsigned int i = 0;
      //cout << line << endl;

      //Extract the name
      string name = "";
      for(; i < line.size(); i++){
        if(line.at(i) == '\t'){
          i++;
          break;
        }
        name += line.at(i);
      }

      //Extract the information flags
      string orientation = "";
      for(; i < line.size(); i++){
        if(line.at(i) == '\t'){
          i++;
          break;
        }
        orientation += line.at(i);
      }

      // Extract the chromosome name
      string chromosome_name = "";
      for(; i < line.size(); i++){
        if(line.at(i) == '\t'){
          i++;
          break;
        }
        chromosome_name += line.at(i);
      }

      // Extract the index in the chromosome
      string index = "";
      for(; i < line.size(); i++){
        if(line.at(i) == '\t'){
          i++;
          break;
        }
        index += line.at(i);
      }


      if(matches->find(name) != matches->end()) {
        unordered_map<string, vector<long>*> * inside_map = matches->at(name);
        if(inside_map->find(chromosome_name) != inside_map->end()){
          inside_map->at(chromosome_name)->push_back(atoi(index.c_str()));
        }
        else{
          vector<long>* v = new vector<long>();
          v->push_back(atoi(index.c_str()));
          inside_map->insert(make_pair(chromosome_name, v));
        }
      }
      else{
        unordered_map<string, vector<long>*> * map = new unordered_map<string, vector<long>*>();
        vector<long> * v = new vector<long>();
        v->push_back(atoi(index.c_str()));
        map->insert(make_pair(chromosome_name, v));
        matches->insert(make_pair(name, map));
      }
    }
    file.close();
  }
  else {
    cerr << "Unable to open file: " << filename << endl;
    exit(0);
  }
}

void read_stdin(unordered_map<string, unordered_map<string, vector<long>* >* > *matches){
  string line;
  while(getline(cin, line)) {
    if(line.size() == 0)
      continue;

    unsigned int i = 0;
    //cout << line << endl;

    //Extract the name
    string name = "";
    for(; i < line.size(); i++){
      if(line.at(i) == '\t'){
        i++;
        break;
      }
      name += line.at(i);
    }

    //cout << name << endl;

    //Extract the information flags
    string orientation = "";
    for(; i < line.size(); i++){
      if(line.at(i) == '\t'){
        i++;
        break;
      }
      orientation += line.at(i);
    }

    // Extract the chromosome name
    string chromosome_name = "";
    for(; i < line.size(); i++){
      if(line.at(i) == '\t'){
        i++;
        break;
      }
      chromosome_name += line.at(i);
    }

    // Extract the index in the chromosome
    string index = "";
    for(; i < line.size(); i++){
      if(line.at(i) == '\t'){
        i++;
        break;
      }
      index += line.at(i);
    }


    if(matches->find(name) != matches->end()) {
      unordered_map<string, vector<long>*> * inside_map = matches->at(name);
      if(inside_map->find(chromosome_name) != inside_map->end()){
        inside_map->at(chromosome_name)->push_back(atoi(index.c_str()));
      }
      else{
        vector<long>* v = new vector<long>();
        v->push_back(atoi(index.c_str()));
        inside_map->insert(make_pair(chromosome_name, v));
      }
    }
    else{
      unordered_map<string, vector<long>*> * map = new unordered_map<string, vector<long>*>();
      vector<long> * v = new vector<long>();
      v->push_back(atoi(index.c_str()));
      map->insert(make_pair(chromosome_name, v));
      matches->insert(make_pair(name, map));
    }
  }
}

void compare_read_matches(unordered_map<string, unordered_map<string, vector<long>* >* >* cxdata, unordered_map<string, unordered_map<string, vector<long>* >* >* aldata){

  // Algorithm: iterate through one map, removing matches in the other map
  // Then iterate through other map, printing the leftover values
  //
  // This algorithm needs to be implemented for the inner and outer maps.

  for(auto i : *cxdata){
    if(aldata->find(i.first) != aldata->end()){
      compare_chromosome_indexes(i.first, i.second, aldata->at(i.first));
      aldata->erase(i.first);
    }
    else{

      cout << i.first << "\t-1\t" << count_false_positives(i.second) << endl;
    }
  }

  for(auto i : *aldata) {
    cout << i.first << "\t1\t0" << endl;
  }

}

void compare_chromosome_indexes(string read_name, unordered_map<string, vector<long>* > * cxdata, unordered_map<string, vector<long>* > * aldata){

  long false_positives = 0;
  string res = "";
  res += read_name;
  res += '\t';
  bool match = false;
  for(auto i : *cxdata) {
    if(aldata->find(i.first) != aldata->end()){
      if(search_vector(i.second, aldata->at(i.first)->at(0))){
        match = true;
        false_positives--;
      }
      false_positives += i.second->size();
    }
    else
      false_positives += i.second->size();
  }
  res += match ? "0\t" : "1\t";
  res += to_string(false_positives);
  cout << res << endl;
}

bool search_vector(vector<long> * source, long x) {
  for(unsigned int i = 0; i < source->size(); i++)
    //if(x==source->at(i))
      //return true;
    //}
    if(abs(x-source->at(i)) <= 5)
      return true;


  return false;
}

long count_false_positives(unordered_map<string, vector<long>* > * cxdata) {
  long false_positives = 0;
  for(auto i : *cxdata) {
    false_positives += i.second->size();
  }
  return false_positives;
}
