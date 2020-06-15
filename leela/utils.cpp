//
//  utils.cpp
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#include "utils.hpp"

//
// functions
//

// find if a suffix exist in string
bool has_suffix(const std::string &str, const std::string &suffix)
{
    return (str.size() >= suffix.size()) && equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

// rename filename if it exist
void filename_rename_ifexist(const std::string &filename)
{
    if (std::ifstream(filename))
    {
        std::string filename_test = filename;
        int counter = 0;
        while (std::ifstream(filename_test))
        {
            filename_test = filename + ".old." + std::to_string(counter);
            ++counter;
        }
        
        if (std::rename(filename.c_str(), filename_test.c_str()))
        {
            std::cerr << filename << " error: " << std::strerror(errno) << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}