#include "Geom.h"

#include <cmath>

static constexpr float PI = 3.14159265359;

namespace geometry
{
    float Distance(float latitude1, float longitude1,
                   float latitude2, float longitude2,
                   bool use_haversine)
    {
        // Haversine is numerically more stable, especially at smaller intervals
        if (!use_haversine)
            return acos(sin(latitude1) * sin(latitude2) + cos(latitude1) * cos(latitude2) * cos(longitude1 - longitude2));

        float del_latitude = latitude2 - latitude2;
        float del_longitude = longitude2 - longitude1;
        float a = sin(del_latitude / 2) * sin(del_latitude / 2) + cos(latitude1) * cos(latitude2) * sin(del_longitude / 2) * sin(del_longitude / 2);
        float c = 2 * atan2(sqrt(a), sqrt(1 - a));
        return c;
    }

    float DegreeMinutesToRealDegrees(int deg, int min)
    {
        constexpr float minute_const = 1.0f / 60;
        float real_degree = (float)deg + minute_const * min;
        return PI * real_degree / 180;
    }

    using namespace data;

    void CalculateDistanceMatrix(const InputDataBuffer& input_data_buffer,
                                 const int num_cities,
                                 Matrix* distance_matrix)
    {
        InitializeMatrix(distance_matrix, num_cities, num_cities);

        for (int i = 0; i < num_cities; i++)
        {
            const Input_Data& city1 = input_data_buffer[i];
            distance_matrix->get(i, i) = 0.0f;

            for (int j = i + 1; j < num_cities; j++)
            {
                const Input_Data& city2 = input_data_buffer[j];

                float dist = Distance(
                    DegreeMinutesToRealDegrees(city1.lat_deg, city1.lat_min),
                    DegreeMinutesToRealDegrees(city1.long_deg, city1.long_min),
                    DegreeMinutesToRealDegrees(city2.lat_deg, city2.lat_min),
                    DegreeMinutesToRealDegrees(city2.long_deg, city2.long_min));

                distance_matrix->get(i, j) = dist;
                distance_matrix->get(j, i) = dist;
            }
        }
    }

    float CalculatePathLength(const Path& path, const Matrix& distance_matrix)
    {
        const int n = path.size();
        float length = 0.0f;

        for (int i = 1; i < n; i++)
        {
            length += distance_matrix.get(path[i], path[i - 1]);
        }

        return length;
    }
} // Namespace geometry