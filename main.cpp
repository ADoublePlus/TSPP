#include "Approx.h"
#include "BranchAndBound.h"
#include "Data.h"
#include "Read.h"

#include <iostream>

using namespace data;
using namespace input;
using namespace geometry;
using namespace branch_and_bound;

template <typename T>
void PrintArr(const std::vector<T>& arr)
{
    for (const T& x : arr)
    {
        std::cout << x << " ";
    }

    std::cout << "\n";
}

static constexpr float EarthRadius = 6371;

void PrintSolution(const System& system, const Path& path)
{
    for (int vert : path)
    {
        std::cout << system.input_data_buffer[vert].city_name << "\n";
    }

    float length = CalculatePathLength(path, system.distance_matrix);

    // NOTE: length is calculated as the degree radians, 
    // multiply it by earth's radius 
    std::cout << length * EarthRadius << " km\n";
}

int main(int argc, char* argv[])
{
    System system;
    system.num_cities = Read("sample.txt", &system.input_data_buffer);
    CalculateDistanceMatrix(system.input_data_buffer, system.num_cities, &system.distance_matrix);
    printf("n: %d\n", system.num_cities);

    // Sydney is the starting city 
    Path path = SolveMetricTSPP<TSPP_Type::FixedStartVert>(system.num_cities, system.distance_matrix, 0, 10000);
    PrintSolution(system, path);
}