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
        ("chr", "Chromosome length file, TSV",   cxxopts::value<string>())
        ("null", "Null model file",               cxxopts::value<string>())
        ("dbsnp",     "SNP database file",             cxxopts::value<string>())
        ("out",       "Output directory",              cxxopts::value<string>()->default_value("./results"))
        ("prefix",    "Output prefix",                 cxxopts::value<string>()->default_value("reli"))
        ("rep",       "Permutation number",            cxxopts::value<int>()->default_value("2000"))
        ("corr",      "Correction multiplier",         cxxopts::value<float>()->default_value("1"));

    opt_parser.parse(argc, argv);

    snps = opt_parser["snp"].as< vector<string> >();
    lds = opt_parser["ld"].as< vector<string> >();
    target = opt_parser["target"].as<string>();
    chr_length = opt_parser["chr"].as<string>();
    null_model = opt_parser["null"].as<string>();
    dbsnp = opt_parser["dbsnp"].as<string>();
    out_dir = opt_parser["out"].as<string>();
    out_prefix = opt_parser["prefix"].as<string>();
    permutation = opt_parser["rep"].as<int>();
    corr_coef = opt_parser["corr"].as<float>();


}

void Parser::print_conf(){
    cout << "--snp && --ld" << endl;
    for (int i = 0; i < snps.size(); i++){
        cout << "  " << snps[i] << endl;
        cout << "  " << lds[i] << endl;
    }
    cout << "--target " << target << endl;
    cout << "--chr " <<chr_length << endl;
    cout << "--null " <<null_model << endl;
    cout << "--dbsnp " <<dbsnp << endl;
    cout << "--out " <<out_dir << endl;
    cout << "--prefix " <<out_prefix << endl;
    cout << "--permutation " <<permutation << endl;
    cout << "--corr " <<corr_coef << endl;
}