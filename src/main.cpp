#include <thread>
#include <vector>
#include <unordered_map>
#include "utils.h"
#include "parser.h"
#include "analysis.h"

using namespace std;


int main(int argc, char* argv[]){

	Parser custom_parser = Parser(argc, argv);

    genomelength genome_length = get_genome_length(custom_parser.chr_length);
    genomelengthsum genome_length_sum = get_genome_length_sum(genome_length);
    snptable snp_table = get_snp_table(custom_parser.dbsnp);

    NullModel null_model(custom_parser.null_model);
    BedData target_bed_data(custom_parser.target);

    vector<statistics> collected_statistics;

    for (int i = 0; i < custom_parser.snps.size(); i++){
        string snp_file = custom_parser.snps[i];
        string ld_file = custom_parser.lds[i];
        cout << "Process: " << endl << "   snp: " << snp_file << endl << "   ld: " << ld_file << endl;
        BedData snp_bed_data(snp_file);
        assign_bins(snp_bed_data, snp_table);
        lddata ld_data = get_ld_data(ld_file, snp_bed_data);
        hits collected_hits = sim(custom_parser.permutation, ld_data, null_model, genome_length, genome_length_sum, target_bed_data);
        statistics stats = get_statistics(collected_hits, ld_data, custom_parser.corr_coef, snp_file);
        collected_statistics.push_back(stats);
    }

    export_results(collected_statistics, custom_parser.out_dir, custom_parser.out_prefix);

    return 0;
}

