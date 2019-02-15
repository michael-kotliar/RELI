/*
Regulatory Element Locus Intersection (RELI) analysis
Copyright (C) <2017>  <Xiaoting.Chen@cchmc.org>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
#include <unistd.h>
#include "RELI_impl.h"
#include "utils.h"
#include "parser.h"
#include "analysis.h"

using namespace std;
using namespace RELI;


int main(int argc, char* argv[]){

	Parser custom_parser = Parser(argc, argv);
    custom_parser.print_conf();
    create_directory(custom_parser.output_dir);

    vector<pair<string, unsigned int> > chr_length_structure = get_chr_length_structure(custom_parser.chrlength_filename);
    vector<unsigned int> chr_length_accumulated = get_chr_length_accumulated(chr_length_structure);

    MafBinnedNullModel null_model_data(custom_parser.nullmodel_filename);
    unordered_map<string, RELI::snp_table_data> snp_table = load_snp_table(custom_parser.dbsnp_filename);

    TargetBedFile target_file(custom_parser.target_filename, null_model_data);

    for (int i = 0; i < custom_parser.snp_filenames.size(); i++){
        RELI::RELIobj *RELIinstance = new RELI::RELIobj;  // dummy object, only for ATGCmap
        string snp_filename = custom_parser.snp_filenames[i];
        string ld_filename = custom_parser.ld_filenames[i];
        vector<SNP> snp_vector = loadSnpFile(snp_filename);
	    extract_snp_info(RELIinstance->ATGCmap, snp_vector, snp_table);
        vector<SNP> snp_vector_temp = snp_vector;
        vector<LD> ld_vector = load_ld_snps(custom_parser.match, ld_filename, snp_vector_temp);
//        vector<double> collected_statistics_vector = sim(null_model_data, ld_vector, custom_parser.match);
//        string output_filename = custom_parser.output_dir + "/" + custom_parser.output_prefix + ".stats";
//        output(output_filename, collected_statistics_vector, ld_vector);
    }
    return 0;
}

