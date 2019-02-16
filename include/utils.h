#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <sys/stat.h>


inline bool create_directory(const string &path){
    cout << "Creating the directory " << path << endl;
    return bool(mkdir(path.c_str(), 0755));
}
#endif