//
// Created by kot4or on 2/14/19.
//

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <unordered_map>
#include <map>
#include <gsl/gsl_cdf.h>
#include <random>
#include <utils.h>


using namespace std;

struct snpline {
    string chr;
    unsigned long start;
    unsigned long end;
    string rsid;
    string obs_strand;
    string ref_allele;
    string alt_alleles;
    string type;
    vector<string> alt_allele_info;
    vector<double> alt_allele_freq;
    snpline(const string &line);
    inline snpline():
            chr(""),
            start(0),
            end(0),
            rsid(""),
            obs_strand(""),
            ref_allele(""),
            alt_alleles(""),
            type("") {};
};


struct atgc{
    map<char, char> mapdata;
    atgc(){
        mapdata['A'] = 'T';
        mapdata['T'] = 'A';
        mapdata['C'] = 'G';
        mapdata['G'] = 'C';
        mapdata['N'] = 'N';
    }
};


struct bed4{
    string chr;
    int uid;
    unsigned long start;
    unsigned long end;
    string name;
    int bin;
    unsigned long length;
    inline bool operator<(const bed4& other) const {
        return (chr < other.chr || (chr == other.chr && start < other.start));
    }
    inline bed4():
            chr(""),
            uid(0),
            start(0),
            end(0),
            name(""),
            bin(0),
            length(0){};
    bed4(const string &line);
    bed4(const bed4 &other_bed4);
};


struct thesamesnp {
    string wanted_snp_id;
    inline bool operator()(const bed4& snp){
        return snp.name == wanted_snp_id;
    }
    thesamesnp(const string& snp_id);
};


struct statistics {
    int intersect;
    unsigned long total;
    double ratio;
    double mean;
    double std;
    double zscore;
    double relative_risk;
    double pval;
    double corrected_pval;
    string source;
    inline statistics():
        intersect(0),
        total(0),
        ratio(0),
        mean(0),
        std(0),
        zscore(0),
        relative_risk(0),
        pval(0),
        corrected_pval(0),
        source(""){}
};



typedef vector<pair<string, unsigned long> > genome;
typedef vector<unsigned long> genomesum;
typedef vector<unsigned long> bin_vector;
typedef unordered_map<string, snpline> snptable;
typedef vector<int> hits;


class NullModel{
public:
    map<int, bin_vector> bin_map;
    NullModel(const string &path);
    void load(const string &path);
};


class BedData{
public:
    vector<bed4> bed_vector;
    map<string, int> bed_index;
    unsigned long median_length;
    BedData(const string &path);
    void load(const string &path);
    void make_index();
};


class LdRecord{
public:
    vector<bed4> bed_vector;  // aka dependent snps
    bed4 key_snp;
    vector<int> distance_to_key_snp;        // Make sure that int is good. And it should be signed
    inline LdRecord(){};
};

typedef vector<LdRecord> lddata;


genome get_genome_structure(const string &path);
genomesum get_genome_sum(const genome &genome_structure);
snptable get_snp_table(const string &path);
lddata get_ld_data(const string &path, const BedData &snp_bed_data);
bool fit_snp (const LdRecord &ld_record, bed4 &temp_snp, const unsigned long &bin, const genome &genome_structure, const genomesum &genome_sum);
hits sim(int permutation, const lddata &ld_data, const NullModel &null_model, const genome &genome_structure, const genomesum &genome_sum, const BedData &target_bed_data);
set<int> get_unique_overlaps(vector<bed4> temp_snp_vector, const BedData &target_bed_data);
bool is_overlap(const bed4 &alpha, const bed4 &beta);
void assign_bins(BedData &snp_bed_data, const snptable &snp_table);
int lookback(int t, int lookback_step = 50);
statistics get_statistics(const hits &collected_hits, const lddata &ld_data, const double &corr_coef,  const string &snp_file, double sig_pct = 0.05);
void export_results(const vector<statistics> &collected_statistics, const string &out_dir, const string &out_prefix);


#endif