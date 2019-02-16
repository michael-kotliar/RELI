//
// Created by kot4or on 2/14/19.

#ifndef PARSER_H
#define PARSER_H

#include <iostream>
#include "cxxopts.hpp"

using namespace std;

class Parser {
public:
    vector<string> snp_filenames;
    vector<string> ld_filenames;
    string target_filename;
    string chrlength_filename;
    string nullmodel_filename;
    string dbsnp_filename;
    string output_dir;
    string output_prefix;
    bool match;
    int permutation;
    int corr_coef;
    string phenotype_name;
    string ancestry_name;
    Parser (int argc, char * argv[]);
    void print_conf();
};

#endif