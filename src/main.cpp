#include <thread>
#include <vector>
#include <map>
#include <unordered_map>
#include "reli.h"
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
        vector<double> collected_statistics_vector = sim(null_model_data, ld_vector, custom_parser.match);
//        string output_filename = custom_parser.output_dir + "/" + custom_parser.output_prefix + ".stats";
//        output(output_filename, collected_statistics_vector, ld_vector);
    }
    return 0;
}

