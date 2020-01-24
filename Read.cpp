#include "Read.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace input
{
    using namespace data;

    int Read(std::string fname, InputDataBuffer* data_buffer)
    {
        std::string line;
        std::string token;
        uint32_t num_cities = 0;
        Input_Data data;

        std::ifstream file(fname);

        while (std::getline(file, line))
        {
            char _comma;
            std::istringstream sline(line);

            // Parse the line and construct the data object
            data.city_id = num_cities;

            std::getline(sline, data.city_name, ',');
            std::getline(sline, token, ',');
            data.lat_deg = atoi(token.c_str());
            std::getline(sline, token, ',');
            data.lat_min = atoi(token.c_str());
            std::getline(sline, token, ',');
            data.long_deg = atoi(token.c_str());
            std::getline(sline, token, ',');
            data.long_min = atoi(token.c_str());

            // Push the data into the buffer 
            data_buffer->push_back(data);

            num_cities++;
        }

        return num_cities;
    }
} // Namespace input