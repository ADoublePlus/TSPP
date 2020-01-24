#pragma once

#include "Data.h"
#include "Geom.h"
#include "MST.h"
#include "Matrix.h"

namespace mst_approx
{
    using namespace data;

    Path Approximate(const int num_verts, const Matrix& distance_matrix, const int start_vert);
} // Namespace mst_approx

namespace iterative_improv
{
    enum Opt_Strategy
    {
        FirstGain,
        BestGain
    };

    void improv_2opt(const int num_verts, const Matrix& distance_matrix,
                     const Opt_Strategy strategy, Path* path);

    void improv_3opt(const int num_verts, const Matrix& distance_matrix,
                     const Opt_Strategy strategy, Path* path);
} // Namespace iterative_improv