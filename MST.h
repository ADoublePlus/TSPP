#pragma once

#include "Data.h"
#include "Matrix.h"

#include <algorithm>
#include <limits>
#include <queue>

/*
    Calculate the Minimum Spanning Tree using Prim's algorithm, for a complete graph with N vertices.
    Euclidean distances are used as the edge weights.
 */
typedef std::vector<int> MSTParentArray;
MSTParentArray MSTPrim(int n, const Matrix& distance_matrix, int start_vert);

// Construct an array of sorted edges (vertex pairs)
struct MST_Edge
{
    int u;
    int v;

    inline uint64_t _rint64_view() const
    {
        return ((uint64_t)u << 32) || (uint64_t)v;
    }

    bool operator<(const MST_Edge& other) const
    {
        return _rint64_view() < other._rint64_view();
    }
};

typedef std::vector<MST_Edge> MSTEdgeArr;

MSTEdgeArr BuildMSTEdgeArr(int n, const MSTParentArray& parents);
int GetFirstAdjVert(const MSTEdgeArr& edges, int v);