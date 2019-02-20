//
// Created by kot4or on 2/16/19.
//

#include <analysis.h>


snpline::snpline(const string &line){
    vector<string> values = split_by_delim(line);
    chr = values[0];
    start = strtoul(values[1].c_str(), NULL, 0);
    end = strtoul(values[2].c_str(), NULL, 0);
    rsid = values[3];
    obs_strand = values[4];
    ref_allele = values[5];
    alt_alleles = values[6];
    type = values[7];
    alt_allele_info = split_by_delim(values[8], ",");
    vector<string> alt_allele_freq_str = split_by_delim(values[9], ",");
    for (int i = 0; i < alt_allele_freq_str.size(); i++)
        alt_allele_freq.push_back(strtod(alt_allele_freq_str[i].c_str(), NULL));
}


NullModel::NullModel(const string &path){
    for (int i = 0; i < 10; i++)
        bin_map[i] = bin_vector();
    load(path);
}

void NullModel::load(const string &path){
    cout << "Loading null model from " << path << endl;
    ifstream input_stream(path);
    string line;
    while (getline(input_stream, line)){
        vector<string> values = split_by_delim(line);
        bin_map[atoi(values[1].c_str())].push_back(strtoul(values[0].c_str(), NULL, 0));
    }
    cout << "- done" << endl;
}

bed4::bed4(const bed4 &other_bed4){
    chr = other_bed4.chr;
    start = other_bed4.start;
    end = other_bed4.end;
    name = other_bed4.name;
    length = other_bed4.length;
    bin = other_bed4.bin;
    uid = other_bed4.uid;
}

bed4::bed4(const string &line){
    vector<string> values = split_by_delim(line);
    chr = values[0];
    start = strtoul(values[1].c_str(), NULL, 0);
    end = strtoul(values[2].c_str(), NULL, 0);
    name = values[3];
    length = end - start;
    bin = 0;
}


BedData::BedData(const string &path){
    load(path);
    make_index();
};


void BedData::make_index(){
    cout << "Calculating index for bed data" << endl;
    string current_chr = bed_vector[0].chr;
    bed_index[current_chr] = 0;
    for (int i = 0; i < bed_vector.size(); i++){
        if (bed_vector[i].chr != current_chr){
            current_chr = bed_vector[i].chr;
            bed_index[current_chr] = i;
        }
    }
    cout << "- done" << endl;
}


void BedData::load(const string &path){
    cout << "Loading bed data from " << path << endl;
    ifstream input_stream(path);
    string line;
    while (getline(input_stream, line)){
        vector<string> values = split_by_delim(line);
        bed4 bed_data(line);
        bed_vector.push_back(bed_data);
    }
    sort(bed_vector.begin(), bed_vector.end());
    vector<unsigned long> length_vector;
    for (int i = 0; i < bed_vector.size(); i++){
        bed_vector[i].uid = i;
        length_vector.push_back(bed_vector[i].length);
    }
    sort(length_vector.begin(), length_vector.end());
    median_length = length_vector[floor(length_vector.size() / 2)];
    cout << "- done" << endl;
}


genome get_genome_structure(const string &path){
    cout << "Loading genome structure from " << path << endl;
    ifstream input_stream(path);
    genome genome_structure;
    string line;
    while(getline(input_stream, line)) {
        vector<string> values = split_by_delim(line);
        pair<string, unsigned long> chr_data (values[0], strtoul(values[1].c_str(), NULL, 0));
        genome_structure.push_back(chr_data);
    }
    cout << "- done" << endl;
    return genome_structure;
}


genomesum get_genome_sum(const genome &genome_structure){
    cout << "Calculating genome structure sum" << endl;
    genomesum genome_sum;
    genome_sum.push_back(0);
    for (int i = 0; i < genome_structure.size(); i++){
        genome_sum.push_back(genome_sum.back() + genome_structure[i].second);
    }
    cout << "- done" << endl;
    return genome_sum;
}


snptable get_snp_table(const string &path){
    cout << "Loading snp table from" << path << endl;
    snptable snp_table;
    ifstream input_stream(path);
    string line;
    while (getline(input_stream, line)){
        snpline snp_line(line);
        snp_table[snp_line.rsid] = snp_line;
    };
    cout << "- done" << endl;
    return snp_table;
}


void assign_bins(BedData &snp_bed_data, const snptable &snp_table){
    cout << "Assigning bins" << endl;
    atgc atgcmap;
    for (auto &bed4_line_it : snp_bed_data.bed_vector){
        auto snp_table_line_iter = snp_table.find(bed4_line_it.name);  // What to do with the snp that are not present in snp_table
        if (snp_table_line_iter != snp_table.end()) {
            vector<double> alt_allele_freq = snp_table_line_iter->second.alt_allele_freq;
            sort(alt_allele_freq.begin(), alt_allele_freq.end());
            bed4_line_it.bin = get_bin(alt_allele_freq.back());
        }
    }
    cout << " - done" << endl;
}


thesamesnp::thesamesnp(const string& snp_id){
    wanted_snp_id = snp_id;
};


lddata get_ld_data(const string &path, const BedData &snp_bed_data){
    cout << "Loading LD data from " << path << endl;
    lddata ld_data;

    ifstream input_stream(path);
    string line;
    map<string, vector<string> > ld_raw_data;
    while (getline(input_stream, line)){
        vector<string> values = split_by_delim(line);
        ld_raw_data[values[0]] = vector<string>(++values.begin(), values.end());
    };

    for (int i = 0; i < snp_bed_data.bed_vector.size(); i++){
        LdRecord ld_record;
        ld_record.key_snp = snp_bed_data.bed_vector[i];
        auto ld_raw_data_iter = ld_raw_data.find(snp_bed_data.bed_vector[i].name);
        if (ld_raw_data_iter != ld_raw_data.end()){
            for (int j = 0; j < ld_raw_data_iter->second.size(); j++){
                auto snp_bed_data_bed_vector_iter = find_if(snp_bed_data.bed_vector.begin(), snp_bed_data.bed_vector.end(), thesamesnp(ld_raw_data_iter->second[j]));
                if (snp_bed_data_bed_vector_iter != snp_bed_data.bed_vector.end()){
                    ld_record.bed_vector.push_back(*snp_bed_data_bed_vector_iter);
                }
            }
        }
        if (ld_record.bed_vector.size() == 0){
            ld_record.bed_vector.push_back(ld_record.key_snp);
        }
        ld_data.push_back(ld_record);
    }

    for (int i = 0; i < ld_data.size(); i++){
        for (int j = 0; j < ld_data[i].bed_vector.size(); j++){
            ld_data[i].distance_to_key_snp.push_back(ld_data[i].bed_vector[j].end - ld_data[i].key_snp.end);
        }
    }
    cout << "- done" << endl;
    return ld_data;
}


bool fit_snp (const LdRecord &ld_record, bed4 &temp_snp, const unsigned long &bin, const genome &genome_structure, const genomesum &genome_sum){

    int max_diff; // Maybe should be unsigned long?
    for (int i = 0; i < ld_record.distance_to_key_snp.size(); i++){
        if (abs(ld_record.distance_to_key_snp[i]) > max_diff){
            max_diff = abs(ld_record.distance_to_key_snp[i]);
        }
    }

    bool found = false;
    for (int i = 0; i < genome_structure.size(); i++){
        if (bin - max_diff - temp_snp.length >= genome_sum[i] && bin + max_diff + temp_snp.length <= genome_sum[i+1]){
            temp_snp.chr = genome_structure[i].first;
            temp_snp.end = bin + (unsigned long)floor(temp_snp.length / 2) - genome_sum[i];
            temp_snp.start = temp_snp.end - temp_snp.length;
            found = true;
            break;
        }
    }

    return found;
}


int lookback(int t, int lookback_step){
    if (t >= lookback_step)
        return (t - lookback_step);
    return 0;
}


bool is_overlap(const bed4 &alpha, const bed4 &beta){
    return (alpha.chr == beta.chr && (
                (alpha.start <= beta.start     && alpha.end       >= beta.end) ||
                (alpha.start >= beta.start     && alpha.start + 1 <  beta.end) ||
                (alpha.end   >  beta.start + 1 && alpha.end       <= beta.end) ||
                (alpha.start >= beta.start     && alpha.end       <= beta.end)));
}


set<int> get_unique_overlaps(vector<bed4> temp_snp_vector, const BedData &target_bed_data){
    set<int> unique_uid_collector;
    int k;
    int t;

    sort(temp_snp_vector.begin(), temp_snp_vector.end());
    string current_chr = temp_snp_vector[0].chr;

    for (int i = 0; i < temp_snp_vector.size(); i++){
        if (temp_snp_vector[i].chr == current_chr){
            t = max(lookback(k), target_bed_data.bed_index.at(temp_snp_vector[i].chr));
            for (k = t; k < target_bed_data.bed_vector.size(); k++){
                if (is_overlap(target_bed_data.bed_vector[k], temp_snp_vector[i])) {
                    unique_uid_collector.insert(temp_snp_vector[i].uid);
                    break;
                }
                if (target_bed_data.bed_vector[k].start >= temp_snp_vector[i].end || target_bed_data.bed_vector[k].chr != temp_snp_vector[i].chr){
                    break;
                }
            }
        } else {
            current_chr = temp_snp_vector[i].chr;
            for (k = target_bed_data.bed_index.at(temp_snp_vector[i].chr); k < target_bed_data.bed_vector.size(); k++){
                if (is_overlap(target_bed_data.bed_vector[k], temp_snp_vector[i])) {
                    unique_uid_collector.insert(temp_snp_vector[i].uid);
                    break;
                }
                if (target_bed_data.bed_vector[k].start >= temp_snp_vector[i].end || target_bed_data.bed_vector[k].chr != temp_snp_vector[i].chr){
                    break;
                }
            }
        }
    }
    return unique_uid_collector;
}


hits sim(int permutation, const lddata &ld_data, const NullModel &null_model, const genome &genome_structure, const genomesum &genome_sum, const BedData &target_bed_data){
    cout << "Running simulation" << endl;
    hits collected_hits;
    default_random_engine random_seed(chrono::system_clock::now().time_since_epoch().count());
    for (int current_iteration = 0; current_iteration < permutation; current_iteration++){
//        cout << "   " << current_iteration + 1 << "/" << permutation << endl;
        lddata simulation_ld_data;
        if (current_iteration == 0){
            simulation_ld_data = ld_data;
        } else {
            for (int current_ld_record = 0; current_ld_record < ld_data.size(); current_ld_record++){
                LdRecord temp_ld_record;
                bed4 temp_snp;
                temp_snp.length = ld_data[current_ld_record].key_snp.length;
                unsigned long t_index;
                bool found = false;
                uniform_int_distribution<unsigned long> dist_generator(0, (null_model.bin_map.at(ld_data[current_ld_record].key_snp.bin).size() - 1));
                while (!found){
                    t_index = dist_generator(random_seed);
                    found = fit_snp(ld_data[current_ld_record], temp_snp, null_model.bin_map.at(ld_data[current_ld_record].key_snp.bin)[t_index], genome_structure, genome_sum);
                }

                for (int current_distance = 0; current_distance < ld_data[current_ld_record].distance_to_key_snp.size(); current_distance++){
                    bed4 second_temp_snp(temp_snp);
                    second_temp_snp.end = second_temp_snp.end + ld_data[current_ld_record].distance_to_key_snp[current_distance];
                    second_temp_snp.start = second_temp_snp.end - second_temp_snp.length;
                    second_temp_snp.uid = current_ld_record;
                    temp_ld_record.bed_vector.push_back(second_temp_snp);
                }
                simulation_ld_data.push_back(temp_ld_record);
            }
        }
        vector<bed4> temp_snp_vector;
        for (int current_simulation_ld_record = 0; current_simulation_ld_record < simulation_ld_data.size(); current_simulation_ld_record++){
            temp_snp_vector.insert(temp_snp_vector.end(), simulation_ld_data[current_simulation_ld_record].bed_vector.begin(), simulation_ld_data[current_simulation_ld_record].bed_vector.end());
        }
        set<int> unique_uid_collector = get_unique_overlaps(temp_snp_vector, target_bed_data);
        collected_hits.push_back(unique_uid_collector.size());
    }
    cout << "- done" << endl;
    return collected_hits;
}


statistics get_statistics(const hits &collected_hits, const lddata &ld_data, const double &corr_coef, const string &snp_file, double sig_pct){
    statistics stats;

    for (int i = 0; i < collected_hits.size(); i++) {
        stats.mean = stats.mean + collected_hits[i];
    }
    stats.mean = stats.mean / double(collected_hits.size());

    double temp = 0;
    for (int i = 0; i < collected_hits.size(); i++) {
        temp += pow((collected_hits[i] - stats.mean), 2);
    }
    stats.std = sqrt(temp / double(collected_hits.size() - 1));

    if (stats.std != 0 && collected_hits[0] >= ceil(double(ld_data.size()) * sig_pct)) {
        stats.zscore = (collected_hits[0] - stats.mean) / stats.std;
    }

    stats.pval = gsl_cdf_ugaussian_Q(stats.zscore);
    stats.corrected_pval = min(stats.pval * corr_coef, 1.0);
    stats.intersect = collected_hits[0];
    stats.total = ld_data.size();
    stats.ratio = stats.intersect / stats.total;
    stats.source = snp_file;
    if (stats.mean != 0){
        stats.relative_risk = stats.intersect / stats.mean;
    }

    return stats;
}

void export_results(const vector<statistics> &collected_statistics, const string &out_dir, const string &out_prefix){
    string path = out_dir + "/" + out_prefix + ".reli";
    cout << "Exporting results to " << path << endl;
    ofstream output_stream (path);
    output_stream // header line
            << "intersect" << "\t"
            << "total" << "\t"
            << "ratio" << "\t"
            << "mean" << "\t"
            << "std" << "\t"
            << "zscore" << "\t"
            << "risk" << "\t"
            << "pval" << "\t"
            << "corr_pval" << "\t"
            << "source" << endl;

    for (int i = 0; i < collected_statistics.size(); i++){
        output_stream
            << collected_statistics[i].intersect << "\t"
            << collected_statistics[i].total << "\t"
            << collected_statistics[i].ratio << "\t"
            << collected_statistics[i].mean << "\t"
            << collected_statistics[i].std << "\t"
            << collected_statistics[i].zscore << "\t"
            << collected_statistics[i].relative_risk << "\t"
            << collected_statistics[i].pval << "\t"
            << collected_statistics[i].corrected_pval << "\t"
            << collected_statistics[i].source << endl;
    }
    output_stream.close();
    cout << " - done" << endl;
}