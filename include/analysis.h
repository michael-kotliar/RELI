//
// Created by kot4or on 2/14/19.
//

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

inline vector<pair<string, unsigned int> > get_chr_length_structure(const string &path){
    ifstream input_stream(path);
    vector<pair<string, unsigned int> > chromosome_structure;
    string line;
    while(getline(input_stream, line)) {
        std::stringstream buffer(line);
        string temp;
        vector<string> values;
        while(getline(buffer, temp, buffer.widen('\t'))) {
            values.push_back(temp);
        }
        pair<string, unsigned int> chr_data (values[0], atoi(values[1].c_str()));
        chromosome_structure.push_back(chr_data);
    }
    return chromosome_structure;
}

inline vector<unsigned int> get_chr_length_accumulated (const vector<pair<string, unsigned int> > &chr_length_structure){
    vector<unsigned int> chr_length_accumulated;
    chr_length_accumulated.push_back(0);
    for (auto k = 0; k < chr_length_structure.size(); ++k){
        chr_length_accumulated.push_back(chr_length_accumulated.back() + chr_length_structure.at(k).second);
    }
    return chr_length_accumulated;
}



#endif
