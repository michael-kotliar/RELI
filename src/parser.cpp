//
// Created by kot4or on 2/14/19.
//

#include <parser.h>

Parser::Parser (int argc, char * argv[]){
    cxxopts::Options opt_parser("RELI", "Regulatory Element Locus Intersection");

    opt_parser.add_options()
        ("snp",       "SNP file",                      cxxopts::value< vector<string> >())
        ("ld",        "LD file",                       cxxopts::value< vector<string> >())
        ("target",    "Target interval file",          cxxopts::value<string>())
        ("chrlength", "Chromosome length file, TSV",   cxxopts::value<string>())
        ("nullmodel", "Null model file",               cxxopts::value<string>())
        ("dbsnp",     "SNP database file",             cxxopts::value<string>())
        ("out",       "Output directory",              cxxopts::value<string>()->default_value("./results"))
        ("prefix",    "Output prefix",                 cxxopts::value<string>()->default_value("reli"))
        ("match",     "Match",                         cxxopts::value<bool>())
        ("rep",       "Permutation number",            cxxopts::value<int>()->default_value("2000"))
        ("corr",      "Correction multiplier",         cxxopts::value<float>()->default_value("1"))
        ("phenotype", "Phenotype name",                cxxopts::value<string>()->default_value("."))
        ("ancestry",  "Ancestry name",                 cxxopts::value<string>()->default_value("."));

    opt_parser.parse(argc, argv);

    snp_filenames = opt_parser["snp"].as< vector<string> >();
    ld_filenames = opt_parser["ld"].as< vector<string> >();
    target_filename = opt_parser["target"].as<string>();
    chrlength_filename = opt_parser["chrlength"].as<string>();
    nullmodel_filename = opt_parser["nullmodel"].as<string>();
    dbsnp_filename = opt_parser["dbsnp"].as<string>();
    output_dir = opt_parser["out"].as<string>();
    output_prefix = opt_parser["prefix"].as<string>();
    permutation = opt_parser["rep"].as<int>();
    corr_coef = opt_parser["corr"].as<float>();
    phenotype_name = opt_parser["phenotype"].as<string>();
    ancestry_name = opt_parser["ancestry"].as<string>();
    match = opt_parser["match"].as<bool>();


}

void Parser::print_conf(){
    cout << "--snp && --ld" << endl;
    for (int i = 0; i < snp_filenames.size(); i++){
        cout << "  " << snp_filenames[i] << endl;
        cout << "  " << ld_filenames[i] << endl;
    }
    cout << "--target " << target_filename << endl;
    cout << "--chrlength " <<chrlength_filename << endl;
    cout << "--nullmodel " <<nullmodel_filename << endl;
    cout << "--dbsnp " <<dbsnp_filename << endl;
    cout << "--out " <<output_dir << endl;
    cout << "--prefix " <<output_prefix << endl;
    cout << "--permutation " <<permutation << endl;
    cout << "--corr " <<corr_coef << endl;
    cout << "--phenotype " <<phenotype_name << endl;
    cout << "--ancestry " <<ancestry_name << endl;
    cout << "--match " <<match << endl;

}