// potential_info.hxx

#pragma once

#include <string>
#include <vector>
#include <iostream>

struct PotentialInfo {
    std::string name;
    std::vector<double> pars;
    int channelid_left;
    int channelid_right;

    // Display potential information
    void summary(std::ostream& os = std::cout) const {
        os << "Potential Information\n";
        os << "---------------------\n";
        os << "Name           : " << name << '\n';
        os << "Channel IDs    : (" << channelid_left << ", " << channelid_right << ")\n";
        os << "Parameters     : [";

        for (std::size_t i = 0; i < pars.size(); ++i) {
            if (i > 0) os << ", ";
            os << pars[i];
        }

        os << "]\n";
    }
};