#pragma once

#include "Matrix.h"

#include <string>
#include <vector>

#define MAX_CITIES 128

namespace data
{
    struct Input_Data
    {
        std::string city_name;
        int city_id = 0;
        int lat_deg = 0;
        int lat_min = 0;
        int long_deg = 0;
        int long_min = 0;
    };

    typedef std::vector<Input_Data> InputDataBuffer;

    struct Output_Data
    {
        std::vector<int> path;
    };

    struct System
    {
        int num_cities;
        InputDataBuffer input_data_buffer;
        Output_Data output_data;
        Matrix distance_matrix;
    };
} // Namespace data

typedef std::vector<int> Path;