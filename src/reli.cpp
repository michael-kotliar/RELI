#include <thread>
#include <mutex>
#include <condition_variable>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <string.h>
#include <time.h>
#include <map>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <unordered_map>
#include <chrono>
#include <cstring>
#include <queue>
#include <random>
#include <assert.h>
#include <unistd.h>
#include "reli.h"
#include "utils.h"
#include <sys/stat.h> // for mkdir and mkdirat (requires kernel >= 2.6.16)
#include <fcntl.h>    // for AT_FDCWD
#include <cerrno>     // for errno


using namespace std;
using namespace RELI;
string buffer;
const int bufferSize = 500000000;
char bufferChar[bufferSize];
static const int lcsize = 500000000;

//bedsig variables
vector<string> RELI::species_name = { "hg19", "mm9" };
vector<RELI::SNP> RELI::SNP_vec;
vector<int> RELI::simulated_number_vec;
vector<RELI::SNP> RELI::SNP_vec_temp;
vector<RELI::LD> RELI::LD_vec;
vector<RELI::LD_template> RELI::LD_template_vec;
vector<RELI::LD_sim> RELI::ldsimvec_after_intersection;
vector<RELI::bed3col> RELI::targetbedinfilevec;
vector<double> RELI::statsvec;
string RELI::TFtype;
string RELI::targetbedinfilename;
string RELI::default_used_species = "hg19";
string RELI::species = RELI::default_used_species;
string RELI::nullmodelinfiledir;
string RELI::defaultBedFileDir;
string RELI::species_chr_mapping_file;
string RELI::nullmodel;
string RELI::bgnullmodel;
string RELI::bgnullmodelinfilename;
string RELI::nullmodelinfilename;
string RELI::dnase_coverage_filename;
string RELI::snpinfilename;
string RELI::temp_nullmodelname, RELI::temp_nullmodelcount;
string RELI::outputpath;
string RELI::outputpath_simulated_intersection;
string RELI::ldfile;
int RELI::secondary_index_size = 50;
RELI::MAF_TSS_binned_null_model RELI::binned_null_model_data2;
unsigned int RELI::dnase_coverage;
RELI::stats_model RELI::used_stats_model = normal;
int RELI::repmax = 2000;
double RELI::corr_muliplier = 1;
double RELI::sig_pct = 0.05;
bool RELI::ldfile_flag = false;
bool RELI::using_default_species = true;
bool RELI::bgnull = false;
double RELI::mu, RELI::sd, RELI::zscore, RELI::pval, RELI::corr_pval;
unordered_map<string, vector<unsigned int>> RELI::indexing_mapping;
map<string, unsigned int> RELI::dnase_coverage_map;
map<string, int> RELI::targetbedfileindex_start;
map<RELI::SNP, RELI::LD_sim, RELI::myless> RELI::snp2ldsim;
RELI::mymap RELI::speciesMap;
vector<pair<string, unsigned int>> RELI::chromosome_strucuture;
vector<unsigned int> RELI::chromosome_strucuture_val;
bool RELI::snp_matching = false;
bool RELI::snp_matching_local_shift = false;
bool RELI::SNPfit(LD LD_A, SNP &tempSNP_A, unsigned int int_A,
	vector<pair<string, unsigned int>> inVec, vector<unsigned int> inVec2){
	int max_diff;
	vector<int> pos_set;
	vector<int> neg_set;
	for (auto it : LD_A.dis2keySNP){
		if (it >0){
			pos_set.push_back(it);
		}
		else{
			neg_set.push_back(it);
		}
	}
	if (pos_set.size()>0){
		max_diff = *max_element(pos_set.begin(), pos_set.end());
	}
	if (neg_set.size() > 0){
		max_diff = max(max_diff, abs(*min_element(neg_set.begin(), neg_set.end())));
	}
	bool okay = false;
	for (auto k = 0; k < inVec.size(); ++k){
		if (int_A - max_diff - tempSNP_A.length >= inVec2.at(k) && int_A + max_diff <= inVec2.at(k + 1)){
			tempSNP_A.snp_chr = inVec.at(k).first;
			tempSNP_A.snp_end = int_A - inVec2.at(k);
			tempSNP_A.snp_start = tempSNP_A.snp_end - tempSNP_A.length;
			okay = true;
			break;
		}
	}
	if (okay == true){
		return true;
	}
	else {
		return false;
	}
}
bool RELI::SNPfit(LD LD_A,
	SNP &tempSNP_A,
	unsigned int int_A,
	vector<pair<string, unsigned int>> inVec,
	vector<unsigned int> inVec2,
	bool inflag){
	int max_diff;
	vector<int> pos_set;
	vector<int> neg_set;
	for (auto it : LD_A.dis2keySNP){
		if (it >0){
			pos_set.push_back(it);
		}
		else{
			neg_set.push_back(it);
		}
	}
	if (pos_set.size()>0){
		max_diff = *max_element(pos_set.begin(), pos_set.end());
	}
	if (neg_set.size() > 0){
		max_diff = max(max_diff, abs(*min_element(neg_set.begin(), neg_set.end())));
	}
	bool okay = false;
	for (auto k = 0; k < inVec.size(); ++k){
		if (int_A - max_diff - tempSNP_A.length >= inVec2.at(k) && int_A + max_diff + tempSNP_A.length <= inVec2.at(k + 1)){
			tempSNP_A.snp_chr = inVec.at(k).first;
			tempSNP_A.snp_end = int_A + floor(tempSNP_A.length / 2) - inVec2.at(k);
			tempSNP_A.snp_start = tempSNP_A.snp_end - tempSNP_A.length;
			okay = true;
			break;
		}
	}
	if (okay == true){
		return true;
	}
	else {
		return false;
	}
}
bool RELI::SNPfit_local(LD &LD_A, //
	SNP &tempSNP_A,		
	unsigned int int_A,
	vector<pair<string, unsigned int>> inVec,
	vector<unsigned int> inVec2,
	bool inflag){
	int max = LD_A.max_dis;
	int min = LD_A.min_dis;
	std::default_random_engine t_randSeed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<int> t_distGen(min, max);
	for (auto &k : LD_A.dis2keySNP){
		k = t_distGen(t_randSeed);
	}
	bool okay = true;
	tempSNP_A.snp_chr = LD_A.keySNP.snp_chr;
	tempSNP_A.snp_start = LD_A.keySNP.snp_end;
	tempSNP_A.snp_end = LD_A.keySNP.snp_end;

	if (okay == true){
		return true;
	}
	else {
		return false;
	}
}
bool RELI::SNPfit_goshift(LD &LD_A,
	SNP &tempSNP_A,		
	unsigned int int_A,
	vector<pair<string, unsigned int>> inVec,
	vector<unsigned int> inVec2,
	bool inflag){
	bool okay = true;
	tempSNP_A.snp_chr = LD_A.keySNP.snp_chr;
	tempSNP_A.snp_start = LD_A.keySNP.snp_end;
	tempSNP_A.snp_end = LD_A.keySNP.snp_end;

	if (okay == true){
		return true;
	}
	else {
		return false;
	}
}
bool RELI::SNPfit_local_index_only(LD &LD_A,
	SNP &tempSNP_A,		
	unsigned int int_A,
	vector<pair<string, unsigned int>> inVec,
	vector<unsigned int> inVec2,
	bool inflag){
	int max = *max_element(LD_A.dis2keySNP.begin(), LD_A.dis2keySNP.end());
	int min = *min_element(LD_A.dis2keySNP.begin(), LD_A.dis2keySNP.end());
	std::default_random_engine t_randSeed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<int> t_distGen(min, max);
	int index_snp_shifting = t_distGen(t_randSeed);
	bool okay = true;
	for (auto i = 0; i < inVec.size(); ++i){
		if (inVec.at(i).first == LD_A.keySNP.snp_chr){
			tempSNP_A.snp_chr = LD_A.keySNP.snp_chr;
			tempSNP_A.snp_start = inVec2.at(i + 1) + LD_A.keySNP.snp_end + index_snp_shifting;
			tempSNP_A.snp_end = inVec2.at(i + 1) + LD_A.keySNP.snp_end + index_snp_shifting;
		}
	}
	if (okay == true){
		return true;
	}
	else {
		return false;
	}
}
void RELI::snpmodifier(SNP &SNP_A, SNP SNP_B, int dist){
	SNP_A.length = SNP_B.length;
	SNP_A.snp_chr = SNP_B.snp_chr;
	SNP_A.snp_end = SNP_B.snp_end + dist;
	SNP_A.snp_start = SNP_A.snp_end - SNP_A.length;
}
int RELI::get_index_to_be_used(string inChr, int inSt, map<pair<string, unsigned int>, int> inMap){
	int rtype;
	vector<pair<unsigned int, int>> tvec;
	for (auto it : inMap){
		if (it.first.first == inChr){
			pair<unsigned int, int> t;
			t.first = it.first.second;
			t.second = it.second;
			tvec.push_back(t);
		}
	}
	if (tvec.size() > 0){
		if (inSt < tvec.begin()->first){
			rtype = tvec.begin()->second;
			return rtype;
		}
		if (inSt >(tvec.end() - 1)->first){
			rtype = (tvec.end() - 1)->second;
			return rtype;
		}
		for (auto it = tvec.begin(); it != tvec.end() - 1; ++it){
			if (inSt >= it->first && inSt < (it + 1)->first){
				rtype = it->second;
				break;
			}
		}
	}
	else{
		rtype = 0;
	}

	return rtype;
}
void RELI::overlapping2(vector<SNP> SNPvecA, vector<bed3col> bedvecA){
	vector<SNP> tempsnpvec;
	tempsnpvec = SNPvecA;
	int k = 0;
	int t;
	sort(tempsnpvec.begin(), tempsnpvec.end(), snpsort);
	string prev_chr = "chr0";

	for (vector<SNP>::iterator snpit = tempsnpvec.begin(); snpit != tempsnpvec.end(); ++snpit){
		if (snpit->snp_chr == prev_chr){
			t = max(lookback_with_zerocheck(k), targetbedfileindex_start[snpit->snp_chr]);
			for (k = t; k < bedvecA.size(); k++){
				if (bedvecA.at(k).bed_chr == snpit->snp_chr &&
                        ((bedvecA.at(k).bed_start <= snpit->snp_start && bedvecA.at(k).bed_end >= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_start + 1 < snpit->snp_end) ||
					(bedvecA.at(k).bed_end > snpit->snp_start + 1 && bedvecA.at(k).bed_end <= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_end <= snpit->snp_end))
					){
					snp2ldsim[*snpit].overlap_sim = true;
					break;
				}
				if (bedvecA.at(k).bed_start >= snpit->snp_end || bedvecA.at(k).bed_chr != snpit->snp_chr){
					break;
				}
			}
		}
		if (snpit->snp_chr != prev_chr){
			prev_chr = snpit->snp_chr;
			for (k = targetbedfileindex_start[snpit->snp_chr]; k <bedvecA.size(); k++){
				if (bedvecA.at(k).bed_chr == snpit->snp_chr &&
                        ((bedvecA.at(k).bed_start <= snpit->snp_start && bedvecA.at(k).bed_end >= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_start + 1 < snpit->snp_end) ||
					(bedvecA.at(k).bed_end > snpit->snp_start + 1 && bedvecA.at(k).bed_end <= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_end <= snpit->snp_end))
					){
					snp2ldsim[*snpit].overlap_sim = true;
					break;
				}
				if (bedvecA.at(k).bed_start >= snpit->snp_end || bedvecA.at(k).bed_chr != snpit->snp_chr){
					break;
				}
			}
		}
	}
}
void RELI::overlapping3(vector<SNP> SNPvecA, vector<bed3col> bedvecA, vector<unsigned int>& in_LD_unique_key_collector){
	vector<SNP> tempsnpvec;
	tempsnpvec = SNPvecA;
	int k = 0;
	int t;
	sort(tempsnpvec.begin(), tempsnpvec.end(), snpsort);
	string prev_chr = "chr0";

	for (vector<SNP>::iterator snpit = tempsnpvec.begin(); snpit != tempsnpvec.end(); ++snpit){
		if (snpit->snp_chr == prev_chr){
			t = max(lookback_with_zerocheck(k), targetbedfileindex_start[snpit->snp_chr]);
			for (k = t; k < bedvecA.size(); k++){
				if (bedvecA.at(k).bed_chr == snpit->snp_chr &&
                        ((bedvecA.at(k).bed_start <= snpit->snp_start && bedvecA.at(k).bed_end >= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_start + 1 < snpit->snp_end) ||
					(bedvecA.at(k).bed_end > snpit->snp_start + 1 && bedvecA.at(k).bed_end <= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_end <= snpit->snp_end))
					){
					in_LD_unique_key_collector.push_back(snpit->inherited_unique_key_from_LD);
					break;
				}
				if (bedvecA.at(k).bed_start >= snpit->snp_end || bedvecA.at(k).bed_chr != snpit->snp_chr){
					break;
				}
			}
		}
		if (snpit->snp_chr != prev_chr){
			prev_chr = snpit->snp_chr;
			for (k = targetbedfileindex_start[snpit->snp_chr]; k <bedvecA.size(); k++){
				if (bedvecA.at(k).bed_chr == snpit->snp_chr &&
                        ((bedvecA.at(k).bed_start <= snpit->snp_start && bedvecA.at(k).bed_end >= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_start + 1 < snpit->snp_end) ||
					(bedvecA.at(k).bed_end > snpit->snp_start + 1 && bedvecA.at(k).bed_end <= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_end <= snpit->snp_end))
					){
					in_LD_unique_key_collector.push_back(snpit->inherited_unique_key_from_LD);
					break;
				}
				if (bedvecA.at(k).bed_start >= snpit->snp_end || bedvecA.at(k).bed_chr != snpit->snp_chr){
					break;
				}
			}
		}
	}
}
void RELI::overlapping_w_index(vector<SNP> SNPvecA, vector<bed3col> bedvecA, vector<unsigned int>& in_LD_unique_key_collector, map<string, int> _index){
	vector<SNP> tempsnpvec;
	tempsnpvec = SNPvecA;
	int k = 0;
	int t;
	sort(tempsnpvec.begin(), tempsnpvec.end(), snpsort);
	string prev_chr = "chr0";

	for (vector<SNP>::iterator snpit = tempsnpvec.begin(); snpit != tempsnpvec.end(); ++snpit){
		if (snpit->snp_chr == prev_chr){
			t = max(lookback_with_zerocheck(k), _index[snpit->snp_chr]);
			for (k = t; k < bedvecA.size(); k++){
				if (bedvecA.at(k).bed_chr == snpit->snp_chr &&
                        ((bedvecA.at(k).bed_start <= snpit->snp_start && bedvecA.at(k).bed_end >= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_start + 1 < snpit->snp_end) ||
					(bedvecA.at(k).bed_end > snpit->snp_start + 1 && bedvecA.at(k).bed_end <= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_end <= snpit->snp_end))
					){
					in_LD_unique_key_collector.push_back(snpit->inherited_unique_key_from_LD);
					break;
				}
				if (bedvecA.at(k).bed_start >= snpit->snp_end || bedvecA.at(k).bed_chr != snpit->snp_chr){
					break;
				}
			}
		}
		if (snpit->snp_chr != prev_chr){
			prev_chr = snpit->snp_chr;
			for (k = _index[snpit->snp_chr]; k <bedvecA.size(); k++){
				if (bedvecA.at(k).bed_chr == snpit->snp_chr &&
                        ((bedvecA.at(k).bed_start <= snpit->snp_start && bedvecA.at(k).bed_end >= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_start + 1 < snpit->snp_end) ||
					(bedvecA.at(k).bed_end > snpit->snp_start + 1 && bedvecA.at(k).bed_end <= snpit->snp_end) ||
					(bedvecA.at(k).bed_start >= snpit->snp_start && bedvecA.at(k).bed_end <= snpit->snp_end))
					){
					in_LD_unique_key_collector.push_back(snpit->inherited_unique_key_from_LD);
					break;
				}
				if (bedvecA.at(k).bed_start >= snpit->snp_end || bedvecA.at(k).bed_chr != snpit->snp_chr){
					break;
				}
			}
		}
	}
}
void RELI::callSpeciesMap(){
	speciesMap["EBNA1"] = "EBV";
	speciesMap["EBNA2"] = "EBV";
	speciesMap["EBNA3C"] = "EBV";
	speciesMap["EBNA3A"] = "EBV";
	speciesMap["EBNA3B"] = "EBV";
	speciesMap["EBNALP"] = "EBV";
	speciesMap["BZLF1"] = "EBV";
	speciesMap["LANA"] = "KSHV";
	speciesMap["TAT"] = "HIV";
}
void RELI::ReadBedSigResult(string inFile, vector<RELI::resultClass>& inVec, bool & inBool){
	ifstream tin;

	tin.open(inFile.c_str());
	while (tin.good()){
		tin.getline(bufferChar, bufferSize);
		tin.peek();
		buffer = bufferChar;
		RELI::resultClass newone;
		newone.track = buffer.substr(0, buffer.find_first_of("\t"));
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.cell = buffer.substr(0, buffer.find_first_of("\t"));
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.tf = buffer.substr(0, buffer.find_first_of("\t"));
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.overlap = atoi(buffer.substr(0, buffer.find_first_of("\t")).c_str());
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.total = atoi(buffer.substr(0, buffer.find_first_of("\t")).c_str());
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.ratio = atof(buffer.substr(0, buffer.find_first_of("\t")).c_str());
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.mean = atof(buffer.substr(0, buffer.find_first_of("\t")).c_str());
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.sd = atof(buffer.substr(0, buffer.find_first_of("\t")).c_str());
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.zscore = atof(buffer.substr(0, buffer.find_first_of("\t")).c_str());
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.enrichment = atof(buffer.substr(0, buffer.find_first_of("\t")).c_str());
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.pval = atof(buffer.substr(0, buffer.find_first_of("\t")).c_str());
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.correctedPval = atof(buffer.substr(0, buffer.find_first_of("\t")).c_str());
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.nullModel = buffer.substr(0, buffer.find_first_of("\t"));
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newone.species = buffer.substr(0, buffer.find_first_of("\t"));

		if (newone.track.find("Rand") != string::npos){
			newone.randomizedSet = true;
			inBool = true;
		}

		inVec.push_back(newone);
	}
	tin.close();
}
void RELI::ReadTFMapping(string inStr, map<string, string>& inMap){
	ifstream tin;

	tin.open(inStr.c_str());
	while (tin.good()){
		tin.getline(bufferChar, bufferSize);
		tin.peek();
		buffer = bufferChar;
		string t = buffer.substr(0, buffer.find_first_of("\t"));
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		inMap[t] = buffer;
	}
	tin.close();
}

void RELI::createSpeciesMap(){
    ifstream tStream;
    tStream.open(RELI::species_chr_mapping_file.c_str());
    while (!tStream.eof()){
        tStream.getline(bufferChar, bufferSize);
        tStream.peek();
        buffer = bufferChar;
        pair<string, unsigned int> t;
        t.first = buffer.substr(0, buffer.find_first_of("\t"));
        buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
        t.second = atoi(buffer.c_str());

        RELI::chromosome_strucuture.push_back(t);
    }
    tStream.close();
	chromosome_strucuture_val.push_back(0);
	for (auto k = 0; k < RELI::chromosome_strucuture.size(); ++k){
		chromosome_strucuture_val.push_back(chromosome_strucuture_val.back() + chromosome_strucuture.at(k).second);
	}
}

bool RELI::mybedsort(bed3col lhs, bed3col rhs){
	return ((lhs.bed_chr < rhs.bed_chr)
		|| ((lhs.bed_chr == rhs.bed_chr) && (lhs.bed_start < rhs.bed_start)));
}

bool RELI::mysort(LD_sim lhs, LD_sim rhs){
	return (lhs.unique_key < rhs.unique_key);
}
bool RELI::myunique(LD_sim lhs, LD_sim rhs){
	return (lhs.unique_key == rhs.unique_key);
}
bool RELI::snpsort(SNP lhs, SNP rhs){
	return ((lhs.snp_chr < rhs.snp_chr) || ((lhs.snp_chr == rhs.snp_chr) && (lhs.snp_start < rhs.snp_start)));
}
unsigned int RELI::mymax(vector<unsigned int> inVec){
	unsigned int maximum = 0;
	for (vector<unsigned int>::iterator intit = inVec.begin(); intit != inVec.end(); ++intit){
		if (*intit >= maximum){ maximum = *intit; };
	}
	return maximum;
}
unsigned int RELI::mymin(vector<unsigned int> inVec){
	unsigned int minimum = 99999999;
	for (vector<unsigned int>::iterator intit = inVec.begin(); intit != inVec.end(); ++intit){
		if (*intit <= minimum){ minimum = *intit; };
	}
	return minimum;
}
int RELI::lookback_with_zerocheck(int t){
	if (t >= RELI::lookback_step){
		return t - lookback_step;
	}
	else{
		return 0;
	}
}
void RELI::cal_stats(RELI::stats_model inModel, const vector<double> &collected_statistics_vector){
	switch (inModel){
	case normal:
	{
		RELI::mu = 0;  // mean
		double temp = 0;
		for (auto it : collected_statistics_vector){
			RELI::mu = RELI::mu + (it);
		}
		RELI::mu = RELI::mu / double(collected_statistics_vector.size());
		for (auto it : collected_statistics_vector){
			temp += (it - RELI::mu)*(it - RELI::mu);
		}
		RELI::sd = sqrt(temp / double(collected_statistics_vector.size() - 1));  // sample sd
		if (RELI::sd == 0 || collected_statistics_vector.at(0) < ceil(double(RELI::LD_vec.size())*RELI::sig_pct)){
			RELI::zscore = 0;
		}
		else{
			RELI::zscore = (collected_statistics_vector.at(0) - RELI::mu) / RELI::sd;
		}
#ifndef bedsig_debug
		RELI::pval = gsl_cdf_ugaussian_Q(RELI::zscore);
#endif
		RELI::corr_pval = min(RELI::pval*RELI::corr_muliplier, 1.0);
	}
	break;

	case empirical:
	{
		auto real_obs = collected_statistics_vector.at(0);
		auto tvec = collected_statistics_vector;
		sort(tvec.begin(), tvec.end());
		int greater_or_equal_instance = 0;
		for (auto i = tvec.rbegin(); i != tvec.rend(); ++i){
			if (*i >= real_obs){
				greater_or_equal_instance++;
			}
		}
		RELI::mu = 0;  // mean
		double temp = 0;
		for (auto it : collected_statistics_vector){
			RELI::mu = RELI::mu + (it);
		}
		RELI::mu = RELI::mu / double(collected_statistics_vector.size());
		for (auto it : collected_statistics_vector){
			temp += (it - RELI::mu)*(it - RELI::mu);
		}
		RELI::sd = sqrt(temp / double(collected_statistics_vector.size() - 1));  // sample sd
		if (RELI::sd == 0 || collected_statistics_vector.at(0) < ceil(double(RELI::LD_vec.size())*RELI::sig_pct)){
			RELI::zscore = 0;
		}
		else{
			RELI::zscore = (collected_statistics_vector.at(0) - RELI::mu) / RELI::sd;
		}
#ifndef bedsig_debug
		RELI::pval = gsl_cdf_ugaussian_Q(RELI::zscore);
#endif
		RELI::corr_pval = double(greater_or_equal_instance) / double(collected_statistics_vector.size());

	}
	break;

	case phasetype:
		break;

	case binomial:
		break;

		// flat the 1000 simulations
	case hypergeometric:

		break;

	case fishers_exact:

		break;

	default:
	{
		RELI::mu = 0;  // mean
		double temp = 0;
		for (auto it : collected_statistics_vector){
			RELI::mu = RELI::mu + (it);
		}
		RELI::mu = RELI::mu / double(collected_statistics_vector.size());
		for (auto it : collected_statistics_vector){
			temp += (it - RELI::mu)*(it - RELI::mu);
		}
		RELI::sd = sqrt(temp / double(collected_statistics_vector.size() - 1));  // sample sd
		if (RELI::sd == 0 || collected_statistics_vector.at(0) < ceil(double(RELI::LD_vec.size())*RELI::sig_pct)){
			RELI::zscore = 0;
		}
		else{
			RELI::zscore = (collected_statistics_vector.at(0) - RELI::mu) / RELI::sd;
		}
#ifndef bedsig_debug
		RELI::pval = gsl_cdf_ugaussian_Q(RELI::zscore);
#endif
		RELI::corr_pval = min(RELI::pval*RELI::corr_muliplier, 1.0);
	}
	break;
	}
}
void RELI::bed3col::cal_avg_peak_length_adjusted_phastCon_score(){
	double tscore = 0;
	double length = this->bed_end - this->bed_start;
	if (length <= 0){
		length = 1;
	}
	for (auto k : this->overlapped_phastCon_data_vec){
		if (k.st > this->bed_start && k.end < this->bed_end){
			tscore += k.sum;
		}
		if (k.st < this->bed_start && k.end < this->bed_end){
			for (auto i = this->bed_start; i < k.end; ++i){
				tscore += k.avg_score;
			}
		}
		if (k.st < this->bed_end && k.end > this->bed_end){
			for (auto i = k.st; i < this->bed_end; ++i){
				tscore += k.avg_score;
			}
		}
		if (k.st<this->bed_start && k.end>this->bed_end){
			if (this->bed_start < this->bed_end){
				for (auto i = this->bed_start; i < this->bed_end; ++i){
					tscore += k.avg_score;
				}
			}
			else{
				tscore += k.avg_score;
			}
		}
	}
	this->avg_peak_length_adjusted_phastCon_score = tscore / length;
}
void RELI::bed3col::cal_avg_peak_length_adjusted_phastCon_score_ez(){
	double bad_value = 0;
	double tscore = 0;
	double length = this->bed_end - this->bed_start;
	if (length <= 0){
		length = 1;
	}
	for (auto k : this->ez_phastCon_score_vec){
		if (k != -1){
			tscore += k;
		}
		else{
			bad_value++;
		}
	}
	if (bad_value < length){
		this->avg_peak_length_adjusted_phastCon_score = tscore / (length - bad_value);
	}
}
void RELI::bed3col::cal_avg_peak_length_adjusted_phastCon_score_ez_50bp(){
	double bad_value = 0;
	double tscore = 0;
	double length = this->ez_phastCon_score_vec_padding_50bp.size();
	if (length <= 0){
		length = 1;
	}
	for (auto k : this->ez_phastCon_score_vec_padding_50bp){
		if (k != -1){
			tscore += k;
		}
		else{
			bad_value++;
		}
	}
	if (bad_value < length){
		this->avg_peak_length_adjusted_phastCon_score = tscore / (length - bad_value);
	}
}
void RELI::TargetBedFile::makeIndex(){
	index[myData.at(0).bed_chr] = 0;
	string prev_chr = myData.at(0).bed_chr;
	for (auto it = myData.begin(); it != myData.end(); ++it){
		if (it->bed_chr != prev_chr){
			prev_chr = it->bed_chr;
			index[it->bed_chr] = distance(myData.begin(), it);
		}
	}
}
void RELI::TargetBedFile::makeIndex2(){
	index[myData.at(0).bed_chr] = 0;
	string prev_chr = myData.at(0).bed_chr;
	for (auto it = myData.begin(); it != myData.end(); ++it){
		if (it->bed_chr != prev_chr){
			int max_increment = (distance(myData.begin(), it) - 1) - index[prev_chr];
			for (auto j = 0; j < secondary_index_size; ++j){
				pair<string, unsigned int> t;
				t.first = it->bed_chr;
				t.second = myData.at(index[prev_chr] + floor(max_increment / secondary_index_size)*j).bed_start;
				this->index2[t] = index[prev_chr] + floor(max_increment / secondary_index_size)*j;
			}
			prev_chr = it->bed_chr;
			index[it->bed_chr] = distance(myData.begin(), it);
		}
	}
}
void RELI::TargetBedFile::cal_avg_phastCons_score(){
	double tscore = 0, tscore_square = 0;
	for (auto k : this->phastCon_score_vec_from_peaks){
		tscore += k;
	}
	this->avg_phastCons_score = tscore / double(this->phastCon_score_vec_from_peaks.size());
	for (auto k : this->phastCon_score_vec_from_peaks){
		tscore_square += (this->avg_phastCons_score - k)*(this->avg_phastCons_score - k);
	}
	this->sd_phastCons_score = sqrt(tscore_square / (double(this->phastCon_score_vec_from_peaks.size() - 1)));
}
void TargetBedFile::load_data(string inStr, bool inVal){
	ifstream in;
	string fname;
	fname = inStr;
	in.open(fname.c_str());
	if (!in){ cout << "can't load file " << fname << endl; exit(1); }
	while (!in.eof()){
		in.getline(bufferChar, bufferSize);
		in.peek();
		buffer = bufferChar;
		RELI::bed3col newTab;
		newTab.bed_chr = buffer.substr(0, buffer.find_first_of("\t"));
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newTab.bed_start = atoi(buffer.substr(0, buffer.find_first_of("\t")).c_str());
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, bufferSize);
		newTab.bed_end = atoi(buffer.substr(0, buffer.find_first_of("\t")).c_str());

		newTab.length = newTab.bed_end - newTab.bed_start;
		this->lengthvec.push_back(newTab.length);
		this->myData.push_back(newTab);
	}
	in.close();
	cout << "target ChIP-seq file loaded, ";
	if (inVal){
		vector<RELI::LD> tLDVec;
		for (auto k : this->myData){
			RELI::LD tld;
			RELI::SNP t;
			t.snp_chr = k.bed_chr;
			t.snp_start = k.bed_start;
			t.snp_end = k.bed_end;
			t.snp_name = "na";
			tld.keySNP = t;
			tld.mySNP.push_back(t);
			tld.dis2keySNP.push_back(0);

			tLDVec.push_back(tld);
		}
		std::default_random_engine tSeed(std::chrono::system_clock::now().time_since_epoch().count()); //RNG seed
		std::uniform_int_distribution<unsigned int> tGen(0, (null_model_data.bin0.size() - 1));// RNG generator
		for (auto k = tLDVec.begin(); k != tLDVec.end(); ++k){
			bool tGood;
			unsigned int tIndex;
			RELI::SNP tKeySNP;
			tKeySNP.length = k->keySNP.length;

			while (tGood != true){
				tIndex = tGen(tSeed);
				tGood = RELI::SNPfit(*k, tKeySNP, null_model_data.bin0.at(tIndex),
					RELI::chromosome_strucuture, RELI::chromosome_strucuture_val);	// fit the new key SNP detail info into simulated location
			}
			bed3col t;
			t.bed_chr = tKeySNP.snp_chr;
			t.bed_start = tKeySNP.snp_start;
			t.bed_end = tKeySNP.snp_end;

			this->myData_bgnull.push_back(t);
			this->myData = this->myData_bgnull;
			this->myData_bgnull.clear();
		}
	}
	sort(this->myData.begin(), this->myData.end());
	cout << "sorted, ";
	sort(this->lengthvec.begin(), this->lengthvec.end());
	this->median_data_length = this->lengthvec.at(floor(this->lengthvec.size() / 2));
	cout << "and width calculated." << endl;
    makeIndex();
}

void RELI::LD::get_features_within_LDblock(const TargetBedFile& rhs){
	map<string, int> local_map_copy = rhs.index;
	for (auto k = local_map_copy[this->LD_chr]; k < rhs.myData.size(); ++k){
		if ((rhs.myData.at(k).bed_chr == this->LD_chr && rhs.myData.at(k).bed_end > this->LD_left_edge && rhs.myData.at(k).bed_end<this->LD_right_edge)
			|| (rhs.myData.at(k).bed_chr == this->LD_chr && rhs.myData.at(k).bed_start> this->LD_left_edge && rhs.myData.at(k).bed_start < this->LD_right_edge)
			|| (rhs.myData.at(k).bed_chr == this->LD_chr && rhs.myData.at(k).bed_start< this->LD_left_edge && rhs.myData.at(k).bed_end > this->LD_right_edge)
			){
			this->features_within_LDblock.myData.push_back(rhs.myData.at(k));
		}
		if (rhs.myData.at(k).bed_chr != this->LD_chr || (rhs.myData.at(k).bed_chr == this->LD_chr && rhs.myData.at(k).bed_start > this->LD_right_edge)){
			break;
		}
	}
}
vector<RELI::bed3col> RELI::LD::goShifting_feature_data(){
	std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<int> dist(0, this->LD_right_edge - this->LD_left_edge);
	vector<RELI::bed3col> tVec;
	for (auto k : this->features_within_LDblock.myData){
		RELI::bed3col tPeak;
		tPeak.bed_chr = k.bed_chr;
		int offset = dist(gen);
		if (offset + k.bed_start > this->LD_right_edge){
			tPeak.bed_start = k.bed_start + offset - this->LD_right_edge + this->LD_left_edge;
			tPeak.bed_end = k.bed_end + offset - this->LD_right_edge + this->LD_left_edge;
		}
		else{
			tPeak.bed_start = k.bed_start + offset;
			tPeak.bed_end = k.bed_end + offset;
		}

		tVec.push_back(tPeak);
	}
	
	return tVec;
}

vector<LD_template> RELI::read_ld_file(const string path){
	vector<LD_template> ld_template_vector;
	ifstream in (path);
	while (!in.eof()){
		in.getline(bufferChar, 5000000);
		buffer = bufferChar;
		buffer.erase(std::remove(buffer.begin(), buffer.end(), '\r'), buffer.end());
		in.peek();

		RELI::LD_template new_ld_template;
		string temp_snp_rs;
		new_ld_template.keySNP = buffer.substr(0, buffer.find_first_of(":"));
		while (buffer.find("\t") != string::npos){
			buffer = buffer.substr(buffer.find_first_of("\t") + 1, 5000000);
			temp_snp_rs = buffer.substr(0, buffer.find_first_of("\t"));
			new_ld_template.mySNP.push_back(temp_snp_rs);
		}
		ld_template_vector.push_back(new_ld_template);   // load into ld template vector
	}
	in.close();
	cout << "phenotype LD file loaded." << endl;
	return ld_template_vector;
}

void RELI::initiate_BedSig(){
	ifstream tin;
	bool match = false;
	for (auto it = RELI::species_name.begin(); it != RELI::species_name.end(); ++it){
		if (*it == RELI::species){
			RELI::nullmodelinfiledir = BedSigNullModelDir + RELI::species + "/";
			RELI::defaultBedFileDir = BedSigBedFileDir + RELI::species + "/";
			RELI::species_chr_mapping_file = BedSigFastaDir + RELI::species + ".bed";
			match = true;
			break;
		}
	}
	if (!match){
		cerr << "no matching species, checking your species with -species option" << endl;
		exit(1);
	}
	//
	RELI::dnase_coverage_filename = RELI::nullmodelinfiledir + "DnaseCoverage.txt";
	tin.open(RELI::dnase_coverage_filename.c_str());
	if (!tin) {
		cerr << "dnase coverage profile not found" << endl;
		exit(1);
	}
	while (!tin.eof()){
		tin.getline(bufferChar, bufferSize);
		tin.peek();
		buffer = bufferChar;

		string buffer0;
		buffer0 = buffer.substr(0, buffer.find_first_of("\t"));
		buffer = buffer.substr(buffer.find_first_of("\t") + 1, buffer.size());
		RELI::dnase_coverage_map[buffer0] = atoi(buffer.c_str());
	}
	tin.close();
	cout << "DnaseCoverage is " << RELI::dnase_coverage << endl;
}
int RELI::binomialCoeff(int n, int k){

	if (k == 0 || k == n){
		return 1;
	}
	else{
		return  binomialCoeff(n - 1, k - 1) + binomialCoeff(n - 1, k);
	}
}
double RELI::binomial_pvalue(int _total, int _overlap, double _prob){
	return double(RELI::binomialCoeff(_total, _overlap))*pow(_prob, double(_overlap)*pow(_prob, double(_total - _overlap)));
}
double RELI::binomial_pvalue_appr(int _total, int _overlap, double _prob){
	return gsl_cdf_ugaussian_Q((_overlap - _total*_prob) / sqrt(_total*_prob*(1 - _prob)));
}
void RELI::RELIobj::public_ver_read_data_index(){
    if (access(this->public_ver_target_label.c_str(), F_OK ) != -1 ) {
        cout << "Skip loading chip-seq index file, because -target points to file" << endl;
    } else {
        ifstream in;
        in.open(this->public_ver_data_index_fname.c_str());
        if (!in){
            cerr << "cannot open data index file, please check with option -index" << endl;
            exit(-1);
        }
        in.ignore(bufferSize, '\n');
        while (!in.eof()){
            in.getline(bufferChar, bufferSize);
            in.peek();
            buffer = bufferChar;
            auto linevec = linehandler(buffer);

            RELI::data_index t;
            t.datalabel = linevec.at(0);
            t.source = linevec.at(1);
            t.cell = linevec.at(2);
            t.tf = linevec.at(3);
            t.cell_label = linevec.at(4);
            t.pmid = linevec.at(5);
            t.group= linevec.at(6);
            t.ebv_status = linevec.at(7);
            t.species = linevec.at(8);

            this->dataindexvec.push_back(t);
        }
        in.close();
        cout << "chip-seq index file loaded." << endl;
    }
}
void RELI::RELIobj::public_ver_set_target_data(){
	if (access(this->public_ver_target_label.c_str(), F_OK ) != -1 ) {
		this->public_ver_target_data_fname = this->public_ver_target_label;
	}
	else {
		auto k = find(this->dataindexvec.begin(),
					  this->dataindexvec.end(),
					  this->public_ver_target_label);
		if (k != this->dataindexvec.end()){
			this->public_ver_selected_data_index = *k;
		}
		else{
			cerr << "cannot find corresponding data entry in the index file, exiting."
				 << endl;
			exit(-1);
		}
		this->public_ver_target_data_fname = this->public_ver_data_dir + "/" + this->public_ver_target_label;
	}
	cout << "target ChIP-seq file set." << endl;
}
void RELI::RELIobj::public_ver_read_null(){
	ifstream in;
}
bool RELI::RELIobj::minimum_check(){
    if (!this->flag_output_prefix){
        this->public_ver_output_prefix = "RELI";
    }

	this->public_ver_output_fname = this->public_ver_output_dir + "/" + this->public_ver_output_prefix + ".stats";
	this->public_ver_output_fname_overlaps = this->public_ver_output_dir + "/" + this->public_ver_output_prefix + ".overlaps";

	cout << "Start Regulatory Element Locus Intersection (RELI) analysis." << endl;
	cout << "Running arguments: " << endl;
	cout << "1) phenotype snp file: " << this->public_ver_snp_fname << endl;								
	cout << "2) phenotype LD structure file: " << RELI::ldfile << endl;								
	cout << "3) SNP matching mode: " << RELI::snp_matching << endl;
	cout << "4) null model file: " << this->public_ver_null_fname << endl;								
	cout << "5) dbSNP table file: " << this->public_ver_snp_table_fname << endl;								
	cout << "6) target chip-seq label / file: " << this->public_ver_target_label << endl;
	cout << "7) chip-seq index file: " << this->public_ver_data_index_fname << endl;								
	cout << "8) chip-seq data dir: " << this->public_ver_data_dir << endl;				
	cout << "9) output dir name: " << this->public_ver_output_dir << endl;
    cout << "10) output file prefix: " << this->public_ver_output_prefix << endl;
	cout << "11) genome build file: " << RELI::species_chr_mapping_file << endl;
	cout << "12) statistic output file name: " << this->public_ver_output_fname << endl;
	cout << "13) overlap output file name: " << this->public_ver_output_fname_overlaps << endl;
	cout << "14) provided phenotype name: " << this->public_ver_phenotype_name << endl;
	cout << "15) provided ancestry name: " << this->public_ver_ancestry_name << endl;
	//
	if (!RELI::snp_matching){
		return (this->flag_input_snp
			&& this->flag_dbsnp_table
			&& this->flag_null_file
			&& this->flag_chipseq_data_dir
			&& this->flag_chipseq_data_index
		
			&& this->flag_output_dir
			&& this->flag_target_label);
	}
	else{
		return (this->flag_input_snp
			&& this->flag_ld_file
			&& this->flag_dbsnp_table
			&& this->flag_null_file
			&& this->flag_chipseq_data_dir
			&& this->flag_chipseq_data_index
			
			&& this->flag_output_dir
			&& this->flag_target_label);
	}
}
void MafBinnedNullModel::load_data(string rhs){
	ifstream in;
	nullmodelinfilename = rhs;
	in.open(nullmodelinfilename.c_str());
	if (!RELI::snp_matching){
		in.ignore(bufferSize, '\n');
	}
	while (!in.eof()){
		in.getline(bufferChar, bufferSize);
		in.peek();
		buffer = bufferChar;
		if (RELI::snp_matching){
			bin_map[atoi(linehandler(buffer).at(1).c_str())]->push_back(atoi(linehandler(buffer).at(0).c_str()));
		}
		else{
			bin0.push_back(atoi(linehandler(buffer).at(0).c_str()));
		}
	}
	in.close();
	cout << "null model loaded." << endl;
}

vector<SNP> RELI::loadSnpFile(const string &path){
	vector<SNP> snp_vector;
	ifstream in (path);
	while (!in.eof()){
		in.getline(bufferChar, 5000);
		in.peek();
		buffer = bufferChar;

		auto buffervec = linehandler(buffer);
		RELI::SNP newsnp;
		newsnp.snp_chr = buffervec.at(0);
		newsnp.snp_start = atoi(buffervec.at(1).c_str());
		newsnp.snp_end = atoi(buffervec.at(2).c_str());
		newsnp.snp_name = buffervec.at(3);

		newsnp.length = newsnp.snp_end - newsnp.snp_start;   // length calculated
		snp_vector.push_back(newsnp);
	}
	in.close();
	return snp_vector;
	cout << "reading snp file completed." << endl;
}

void RELI::RELIobj::create_output_dir(){
	int rv;
	struct stat st;
	const char *outpath_c = this->public_ver_output_dir.c_str();

	// don't try to make the directory if it already exists
	if (stat(outpath_c, &st) < 0) {
		if (mkdirat(AT_FDCWD, outpath_c, 0755) < 0 && errno != EEXIST) {
			cerr << "Unable to create output dir '" << outpath_c << "': "
			     << strerror(errno) << "." << endl;
			exit(1);
		}
	} else if (!S_ISDIR(st.st_mode)) {
		cerr << "Specified output dir '" << outpath_c
		     << "' is not a directory!" << endl;
		exit(1);
	}
}

unordered_map<string, RELI::snp_table_data> RELI::load_snp_table(const string &path){
    unordered_map<string, RELI::snp_table_data> snp_table;
	ifstream in(path);
	in.ignore(lcsize,'\n');
	while (!in.eof()){
		in.getline(bufferChar, bufferSize);
		in.peek();
		buffer = bufferChar;
		auto linevec = linehandler(buffer);

		RELI::snp_table_data t;
		t.chr = linevec.at(0);
		t.start = linevec.at(1);
		t.end = linevec.at(2);
		t.rsid = linevec.at(3);
		t.obs_strand = linevec.at(4);
		t.ref_allele = linevec.at(5);
		t.alt_alleles = linevec.at(6);
		t.type= linevec.at(7);
		t.alt_allele_info = linevec.at(8);
		t.alt_allele_freq = linevec.at(9);
        snp_table[t.rsid] = t;
	};
	in.close();
    return snp_table;
}

string RELI::dnaSeqReverse(string inSeq, map<char, char> thismap){
	string tempseq;
	for (string::reverse_iterator sit = inSeq.rbegin(); sit != inSeq.rend(); ++sit){
		tempseq += thismap[*sit];
	}
	return tempseq;
}

void RELI::extract_snp_info(map<char,char> rhs, vector<SNP> &snp_vector, unordered_map<string, RELI::snp_table_data> &snp_table){
	for (auto &snp_it : snp_vector){
		if (snp_table[snp_it.snp_name].chr.size()>0){
			snp_it.obs_strand = snp_table[snp_it.snp_name].obs_strand;
			snp_it._ref_allele = snp_table[snp_it.snp_name].ref_allele;
			string alt_alleles = snp_table[snp_it.snp_name].alt_alleles;
			snp_it.snp_type = snp_table[snp_it.snp_name].type;
			if (snp_it.obs_strand == "-"){
				snp_it._alt_allele.push_back(RELI::dnaSeqReverse(alt_alleles.substr(0, alt_alleles.find_first_of("/")), rhs));
				while (alt_alleles.find("/") != string::npos){
					alt_alleles = alt_alleles.substr(alt_alleles.find_first_of("/") + 1, alt_alleles.size());
					snp_it._alt_allele.push_back(RELI::dnaSeqReverse(alt_alleles.substr(0, alt_alleles.find_first_of("/")), rhs));
				}
			}
			else{
				snp_it._alt_allele.push_back(alt_alleles.substr(0, alt_alleles.find_first_of("/")));
				while (alt_alleles.find("/") != string::npos){
					alt_alleles = alt_alleles.substr(alt_alleles.find_first_of("/") + 1, alt_alleles.size());
					snp_it._alt_allele.push_back(alt_alleles.substr(0, alt_alleles.find_first_of("/")));
				}
			}
			string alt_alleles_2 = snp_table[snp_it.snp_name].alt_allele_info;
			string alt_alleles_freq = snp_table[snp_it.snp_name].alt_allele_freq;
			snp_it._MAF_alt_allele_string = alt_alleles_2;
			if (alt_alleles_2.size() > 0 && alt_alleles_freq.size()>0){
				vector<double> tvec;
				while (alt_alleles_freq.find(",") != string::npos){
					tvec.push_back(atof(alt_alleles_freq.substr(0, alt_alleles_freq.find_first_of(",")).c_str()));
					alt_alleles_freq = alt_alleles_freq.substr(alt_alleles_freq.find_first_of(",") + 1, alt_alleles_freq.size());
				}
				sort(tvec.begin(), tvec.end());
				snp_it._MAF = *(tvec.end() - 2);
				snp_it._MAF_Bin = snp_it.cal_MAF_Bin(snp_it._MAF);
			}
			else{
				snp_it._MAF = 0.001;
				snp_it._MAF_Bin = snp_it.cal_MAF_Bin(snp_it._MAF);
			}
		}
		else{
			snp_it._MAF = 0.001;
			snp_it._MAF_Bin = snp_it.cal_MAF_Bin(snp_it._MAF);
		}
	}
	cout << "snp MAF information queried." << endl;
}
vector<LD> RELI::load_ld_snps(bool match, string ld_path, vector<SNP> &snp_vector_temp){
	vector<LD> ld_vector;
	ifstream in;
	if (match){
		vector<LD_template> ld_template_vector = RELI::read_ld_file(ld_path);
		for (auto it = ld_template_vector.begin(); it != ld_template_vector.end(); ++it){
			RELI::LD newld;  // real LD instance
			newld.keySNP = *find(snp_vector_temp.begin(), snp_vector_temp.end(), it->keySNP);
			for (auto k = 0; k < it->mySNP.size(); k++){
				for (auto snpit = snp_vector_temp.begin(); snpit != snp_vector_temp.end(); ++snpit){
					if (it->mySNP.at(k) == snpit->snp_name){
						newld.mySNP.push_back(*snpit);
						snp_vector_temp.erase(snpit);
						break;
					}
				}
			}
			ld_vector.push_back(newld);
		}
		for (auto snpit = snp_vector_temp.begin(); snpit != snp_vector_temp.end(); ++snpit){
			RELI::LD newld;
			newld.keySNP = *snpit;
			newld.mySNP.push_back(*snpit);

			ld_vector.push_back(newld);
		}
	}
	else{
		for (auto snpit = snp_vector_temp.begin(); snpit != snp_vector_temp.end(); ++snpit){
			RELI::LD newld;
			newld.keySNP = *snpit;
			newld.mySNP.push_back(*snpit);

			ld_vector.push_back(newld);
		}
	}
	for (auto ldit = ld_vector.begin(); ldit != ld_vector.end(); ++ldit){
		for (auto snpit = ldit->mySNP.begin(); snpit != ldit->mySNP.end(); ++snpit){
			ldit->dis2keySNP.push_back(snpit->snp_end - ldit->keySNP.snp_end);
		}
	}
	cout << "LD structure handled. " << endl;
	return ld_vector;
}

void RELI::output(const string &path, const vector<double> &collected_statistics_vector, const vector<LD> &ld_vector){
	RELI::cal_stats(RELI::used_stats_model, collected_statistics_vector);
	ofstream out(path);
	out << "Formal Phenotype" << "\t" << "Ancestry" << "\t" << "Source"
		<< "\t" << "Cell" << "\t" << "Formal Cell" << "\t" << "Label"
		<< "\t" << "Intersect" << "\t" << "Total" << "\t" << "Ratio"
		<< "\t" << "Mean" << "\t" << "Std" << "\t" << "Z-score" << "\t"
		<< "Relative Risk" << "\t" << "P-val" << "\t" << "Corrected P-val"
		<< "\t" << "Null_Model" << "\t" << "Species" << endl;
	out << "\t" <<  collected_statistics_vector[0]
		<< "\t" << ld_vector.size()
		<< "\t" << collected_statistics_vector[0] / ld_vector.size()
		<< "\t" << RELI::mu
		<< "\t" << RELI::sd
		<< "\t" << RELI::zscore
		<< "\t";
		if (RELI::mu != 0){
			out << collected_statistics_vector[0] / RELI::mu;
		}
		else{
			out << 0;
		}
	out << "\t" << RELI::pval
		<< "\t" << RELI::corr_pval
		<< endl;
	out.close();
	cout << "RELI analysis completed, check output folder for exciting discoveries!" << endl;
}

vector<double> RELI::sim(MafBinnedNullModel null_model, const vector<LD> &ld_vector, bool match){
	std::default_random_engine randSeed(std::chrono::system_clock::now().time_since_epoch().count());
	vector<double> collected_statistics_vector;
	for (auto i = 0; i < repmax + 1; i++){
		vector<RELI::LD_sim> ld_simulation_vector;
		if (i == 0){
			for (auto LDit = ld_vector.begin(); LDit != ld_vector.end(); ++LDit){
				RELI::LD_sim t_LD_sim;
				t_LD_sim.unique_key = distance(ld_vector.begin(), LDit);
				t_LD_sim.mySNP = LDit->mySNP;   	
				for (auto &tSNP_iter : t_LD_sim.mySNP){
					tSNP_iter.inherited_unique_key_from_LD = t_LD_sim.unique_key;
				}
				t_LD_sim.dis2keySNP = LDit->dis2keySNP;
				t_LD_sim.keySNP = LDit->keySNP;

				t_LD_sim.overlap_sim = false;
				ld_simulation_vector.push_back(t_LD_sim);
			}
		}
		if (i>0){
			for (auto LDit = ld_vector.begin(); LDit != ld_vector.end(); ++LDit){
				RELI::LD_sim t_LD_sim;
				t_LD_sim.unique_key = distance(ld_vector.begin(), LDit);
				unsigned int tIndex;
				RELI::SNP tKeySNP;
				tKeySNP.length = LDit->keySNP.length;
				bool datagood = false;
				if (match){
					std::uniform_int_distribution<unsigned int> distGen(0, (null_model.bin_map[LDit->keySNP._MAF_Bin]->size() - 1));
					while (datagood != true){
						tIndex = distGen(randSeed);
						datagood = RELI::SNPfit(*LDit,
							tKeySNP,
                            null_model.bin_map[LDit->keySNP._MAF_Bin]->at(tIndex),
							RELI::chromosome_strucuture,
							RELI::chromosome_strucuture_val,
							1);	
					}
				}
				else{
					std::uniform_int_distribution<unsigned int> distGen(0, (null_model.bin0.size() - 1));
					while (datagood != true){
						tIndex = distGen(randSeed);
						datagood = RELI::SNPfit(*LDit,
							tKeySNP, null_model.bin0.at(tIndex),
							RELI::chromosome_strucuture,
							RELI::chromosome_strucuture_val,
							1);	
					}
				}
				for (auto Dit = LDit->dis2keySNP.begin(); Dit != LDit->dis2keySNP.end(); ++Dit){
					RELI::SNP tSNP;
					snpmodifier(tSNP, tKeySNP, *Dit);
					tSNP.inherited_unique_key_from_LD = t_LD_sim.unique_key;
					t_LD_sim.mySNP.push_back(tSNP);
				}
				t_LD_sim.overlap_sim = false;
				ld_simulation_vector.push_back(t_LD_sim);
			}
		}
		vector<SNP> snp_vector_temp;
		for (auto ldsimit = ld_simulation_vector.begin(); ldsimit != ld_simulation_vector.end(); ++ldsimit){
			for (auto snpit = ldsimit->mySNP.begin(); snpit != ldsimit->mySNP.end(); ++snpit){
				snp_vector_temp.push_back(*snpit);
			}
		}
		vector<unsigned int> LD_unique_key_collector;
		RELI::overlapping3(snp_vector_temp, RELI::targetbedinfilevec, LD_unique_key_collector);   // only updated indicator in mapped ldsim instances
		sort(LD_unique_key_collector.begin(), LD_unique_key_collector.end());
		LD_unique_key_collector.resize(distance(LD_unique_key_collector.begin(), unique(LD_unique_key_collector.begin(), LD_unique_key_collector.end())));
		collected_statistics_vector.push_back(double(LD_unique_key_collector.size()));
		if (i % 500 == 0){
			cout << double(i)/double(RELI::repmax)*100<<"% finished." << endl;
		}
	}
	return collected_statistics_vector;
}

vector<string> RELI::linehandler(string inStr){
	vector<string> robj;
	robj.push_back(inStr.substr(0, inStr.find_first_of("\t")));
	while (inStr.find("\t") != string::npos){
		inStr = inStr.substr(inStr.find_first_of("\t") + 1, inStr.size());
		robj.push_back(inStr.substr(0, inStr.find_first_of("\t")));
	}
	return robj;
}
