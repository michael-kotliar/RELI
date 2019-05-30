//
// Created by kot4or on 2/16/19.
//

#include <analysis.h>


snprecord::snprecord(const string &line){
    vector<string> values = split_by_delim(line);
    chr = values[0];
    start = strtoul(values[1].c_str(), NULL, 0);
    end = strtoul(values[2].c_str(), NULL, 0);
    rsid = values[3];
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
    cout << endl <<"Loading null model from " << path << endl;
    ifstream input_stream(path);
    string line;
    while (getline(input_stream, line)){
        vector<string> values = split_by_delim(line);
        bin_map[atoi(values[1].c_str())].push_back(strtoul(values[0].c_str(), NULL, 0));
    }
    for (auto it = bin_map.begin(); it != bin_map.end(); ++it){
        cout << it->first << ": " << it->second.size() << endl;
    }
}

bed4::bed4(const bed4 &other_bed4){
    chr = other_bed4.chr;
    name = other_bed4.name;
    start = other_bed4.start;
    end = other_bed4.end;
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
    cout << endl << "Calculating index for bed data" << endl;
    string current_chr = bed_vector[0].chr;
    bed_index[current_chr] = 0;
    for (int i = 0; i < bed_vector.size(); i++){
        if (bed_vector[i].chr != current_chr){
            current_chr = bed_vector[i].chr;
            bed_index[current_chr] = i;
        }
    }
    
    for (auto &bed_index_it : bed_index){
        cout << setw(5) << bed_index_it.first << ": " << setw(10) << bed_index_it.second << endl;
    }
}


void BedData::load(const string &path){
    cout << endl << "Loading bed data from " << path << endl;
    ifstream input_stream(path);
    string line;
    while (getline(input_stream, line)){
        bed4 bed_data(line);
        bed_vector.push_back(bed_data);
    }

    cout << endl << "Sorting" << endl << endl;
    sort(bed_vector.begin(), bed_vector.end());

    set<string> uniq_chr_set;
    for (int i = 0; i < bed_vector.size(); i++){
        if (uniq_chr_set.find(bed_vector[i].chr) == uniq_chr_set.end()){
            cout << bed_vector[i].chr << endl;
            uniq_chr_set.insert(bed_vector[i].chr);
        }
    }

    vector<unsigned long> length_vector;
    for (int i = 0; i < bed_vector.size(); i++){
        bed_vector[i].uid = i;
        length_vector.push_back(bed_vector[i].length);
    }
    sort(length_vector.begin(), length_vector.end());
    median_length = length_vector[floor(length_vector.size() / 2)];
}


genomelength get_genome_length(const string &path){
    cout << endl << "Loading genome length from " << path << endl << endl;
    ifstream input_stream(path);
    genomelength genome_length;
    string line;
    while(getline(input_stream, line)) {
        vector<string> values = split_by_delim(line);
        pair<string, unsigned long> chr_data (values[0], strtoul(values[1].c_str(), NULL, 0));
        cout << setw(5) << chr_data.first << ": " << setw(10) << chr_data.second << endl;
        genome_length.push_back(chr_data);
    }
    
    cout << endl << "Sorting" << endl << endl;
    sort(genome_length.begin(), genome_length.end());
    
    for (int i = 0; i < genome_length.size(); i++){
        cout << setw(5) << genome_length[i].first << ": " << setw(10) << genome_length[i].second << endl;
    }

    return genome_length;
}


genomelengthsum get_genome_length_sum(const genomelength &genome_structure){
    cout << endl << "Calculating genome length sum" << endl << endl;
    genomelengthsum genome_length_sum;
    genome_length_sum.push_back(0);
    cout << setw(27) << " | " << genome_length_sum.back() << endl;
    for (int i = 0; i < genome_structure.size(); i++){
        genome_length_sum.push_back(genome_length_sum.back() + genome_structure[i].second);
        cout << setw(5) << genome_structure[i].first << ": " << setw(10) << genome_structure[i].second << setw(10) << " | " << genome_length_sum.back() << endl;
    }
    return genome_length_sum;
}


snptable get_snp_table(const string &path){
    cout << endl << "Loading snp table from" << path << endl << endl;
    snptable snp_table;
    ifstream input_stream(path);
    string line;
    int i = 0;
    while (getline(input_stream, line)){
        snprecord snp_record(line);
        snp_table[snp_record.rsid] = snp_record;
        if (i < 10){
            cout << setw(15) << snp_record.rsid
                 << ": "  << setw(10) << snp_table[snp_record.rsid].chr
                 << " "   << setw(10) << snp_table[snp_record.rsid].start
                 << " "   << setw(10) << snp_table[snp_record.rsid].end
                 << " | " << setw(10);
            for (int j = 0; j < snp_table[snp_record.rsid].alt_allele_freq.size(); j++){
                cout << " " << setw(10) << snp_table[snp_record.rsid].alt_allele_freq[j];
            }
            cout << endl;
            i++;
        }
    };
    return snp_table;
}


void assign_bins(BedData &snp_bed_data, const snptable &snp_table){
    cout << endl << "Assigning bins" << endl << endl;
    atgc atgcmap;
    for (auto &bed4_line_it : snp_bed_data.bed_vector){
        auto snp_table_line_iter = snp_table.find(bed4_line_it.name);  // What to do with the snp that are not present in snp_table
        if (snp_table_line_iter != snp_table.end()) {
            vector<double> alt_allele_freq = snp_table_line_iter->second.alt_allele_freq;
            sort(alt_allele_freq.begin(), alt_allele_freq.end());
            bed4_line_it.bin = get_bin(alt_allele_freq[alt_allele_freq.size() - 2]); // second from the last item
        }
    }
}


thesamesnp::thesamesnp(const string& snp_id){
    wanted_snp_id = snp_id;
};


lddata get_ld_data(const string &path, const BedData &snp_bed_data){
    cout << endl << "Loading LD data from " << path << endl << endl;
    lddata ld_data;

    ifstream input_stream(path);
    string line;
    map<string, vector<string> > ld_raw_data;
    while (getline(input_stream, line)){
        vector<string> values = split_by_delim(line);
        ld_raw_data[values[0]] = vector<string>(++values.begin(), values.end());
    };


    for (auto &ld_raw_data_it : ld_raw_data){
        cout << setw(5) << ld_raw_data_it.first << ": ";
        for (int i = 0; i <  ld_raw_data_it.second.size(); i++){
            cout << setw(20) << ld_raw_data_it.second[i];
        }
        cout << endl;

        LdRecord ld_record;
        ld_record.key_snp = *find_if(snp_bed_data.bed_vector.begin(), snp_bed_data.bed_vector.end(), thesamesnp(ld_raw_data_it.first));
        for (int j = 0; j < ld_raw_data_it.second.size(); j++){
            auto snp_bed_data_bed_vector_iter = find_if(snp_bed_data.bed_vector.begin(), snp_bed_data.bed_vector.end(), thesamesnp(ld_raw_data_it.second[j]));
            if (snp_bed_data_bed_vector_iter != snp_bed_data.bed_vector.end()){
                ld_record.bed_vector.push_back(*snp_bed_data_bed_vector_iter);
            } else {
                ld_record.bed_vector.push_back(ld_record.key_snp);  // looks like bypass. It should never happen if LD was generated based on SNP properly
            }
        }
        ld_data.push_back(ld_record);
    }

    // for (int i = 0; i < snp_bed_data.bed_vector.size(); i++){
    //     LdRecord ld_record;
    //     ld_record.key_snp = snp_bed_data.bed_vector[i];
    //     auto ld_raw_data_iter = ld_raw_data.find(ld_record.key_snp.name);
    //     if (ld_raw_data_iter != ld_raw_data.end()){
    //         for (int j = 0; j < ld_raw_data_iter->second.size(); j++){
    //             auto snp_bed_data_bed_vector_iter = find_if(snp_bed_data.bed_vector.begin(), snp_bed_data.bed_vector.end(), thesamesnp(ld_raw_data_iter->second[j]));
    //             if (snp_bed_data_bed_vector_iter != snp_bed_data.bed_vector.end()){
    //                 ld_record.bed_vector.push_back(*snp_bed_data_bed_vector_iter);
    //             }
    //         }
    //     }
    //     if (ld_record.bed_vector.size() == 0){
    //         ld_record.bed_vector.push_back(ld_record.key_snp);
    //     }
    //     ld_data.push_back(ld_record);
    // }

    for (int i = 0; i < ld_data.size(); i++){
        for (int j = 0; j < ld_data[i].bed_vector.size(); j++){
            ld_data[i].distance_to_key_snp.push_back(ld_data[i].bed_vector[j].end - ld_data[i].key_snp.end);
        }
    }
    
    cout << endl << "Calculated distances to SNPs from SNP file" << endl << endl;
    for (int i = 0; i < ld_data.size(); i++){
        cout << setw(5) << ld_data[i].key_snp.name << ": ";
        for (int j = 0; j <  ld_data[i].distance_to_key_snp.size(); j++){
            cout << setw(10) << ld_data[i].distance_to_key_snp[j];
        }
        cout << endl;
    }

    return ld_data;
}


bool fit_snp (const LdRecord &ld_record, bed4 &temp_snp, const unsigned long &bin, const genomelength &genome_structure, const genomelengthsum &genome_sum){

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


hits sim(int permutation, const lddata &ld_data, const NullModel &null_model, const genomelength &genome_structure, const genomelengthsum &genome_sum, const BedData &target_bed_data){
    cout << endl << "Running simulation" << endl << endl;
    hits collected_hits;
    default_random_engine random_seed(chrono::system_clock::now().time_since_epoch().count());
    for (int current_iteration = 0; current_iteration < permutation; current_iteration++){
        cout << "   " << current_iteration + 1 << "/" << permutation << endl;
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
    return collected_hits;
}


statistics get_statistics(const hits &collected_hits, const lddata &ld_data, const double &corr_coef, const string &snp_file, double sig_pct){
    statistics stats;

    for (int i = 0; i < collected_hits.size(); i++) {
        stats.mean = double(stats.mean) + double(collected_hits[i]);
    }
    stats.mean = double(stats.mean) / double(collected_hits.size());

    double temp = 0;
    for (int i = 0; i < collected_hits.size(); i++) {
        temp += pow((double(collected_hits[i]) - double(stats.mean)), 2);
    }
    stats.std = sqrt(double(temp) / double(collected_hits.size() - 1));

    if (stats.std != 0 && collected_hits[0] >= ceil(double(ld_data.size()) * double(sig_pct))) {
        stats.zscore = (double(collected_hits[0]) - double(stats.mean)) / double(stats.std);
    }

    stats.pval = double(gsl_cdf_ugaussian_Q(double(stats.zscore)));
    stats.corrected_pval = min(double(stats.pval) * double(corr_coef), 1.0);
    stats.intersect = double(collected_hits[0]);
    stats.total = double(ld_data.size());
    stats.ratio = double(stats.intersect) / double(stats.total);
    stats.source = snp_file;
    if (stats.mean != 0){
        stats.relative_risk = double(stats.intersect) / double(stats.mean);
    }

    return stats;
}

void export_results(const vector<statistics> &collected_statistics, const string &out_dir, const string &out_prefix){
    create_directory(out_dir);
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
}