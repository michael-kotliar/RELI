#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <sys/stat.h>
#include <sstream>
#include <iostream>

using namespace std;


inline bool create_directory(const string &path){
    cout << "Creating output directory " << path << endl;
    cout << "- done" << endl;
    return bool(mkdir(path.c_str(), 0755));
}


inline vector<string> split_by_delim(const string &line, string delim = "\t"){
    stringstream buffer(line);
    string temp;
    vector<string> values;
    while(getline(buffer, temp, buffer.widen(*delim.c_str()))) {
        values.push_back(temp);
    }
    return values;
}


inline int get_bin(const double &maf){
    if (0.05 <= maf && maf < 0.1)  { return 1; }
    if (0.1  <= maf && maf < 0.15) { return 2; }
    if (0.15 <= maf && maf < 0.2)  { return 3; }
    if (0.2  <= maf && maf < 0.25) { return 4; }
    if (0.25 <= maf && maf < 0.3)  { return 5; }
    if (0.3  <= maf && maf < 0.35) { return 6; }
    if (0.35 <= maf && maf < 0.4)  { return 7; }
    if (0.4  <= maf && maf < 0.45) { return 8; }
    if (0.45 <= maf && maf <= 0.5) { return 9; }
    return 0;
}


#endif