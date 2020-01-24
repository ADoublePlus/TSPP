#pragma once

#include "Data.h"
#include "Matrix.h"

namespace geometry
{
    float Distance(float latitude1, float longitude1,
                   float latitude2, float longitude2,
                   bool use_haversine = true);

    float DegreeMinutesToRealDegrees(uint8_t deg, uint8_t min);

    using namespace data;

    void CalculateDistanceMatrix(const InputDataBuffer& input_data_buffer,
                                 const int num_cities,
                                 Matrix* distance_matrix);

    float CalculatePathLength(const Path& path, const Matrix& distance_matrix);
} // Namespace geometry