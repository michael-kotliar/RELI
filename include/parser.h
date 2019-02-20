//
// Created by kot4or on 2/14/19.

#ifndef PARSER_H
#define PARSER_H

#include <iostream>
#include "cxxopts.hpp"

using namespace std;

class Parser {
public:
    vector<string> snps;
    vector<string> lds;
    string target;
    string chr_length;
    string null_model;
    string dbsnp;
    string out_dir;
    string out_prefix;
    int permutation;
    int corr_coef;
    string phenotype_name;
    string ancestry_name;
    Parser (int argc, char * argv[]);
    void print_conf();
};

#endif