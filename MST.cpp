#include "MST.h"
#include "Flag_Array.h"

MSTParentArray MSTPrim(int n, const Matrix& distance_matrix, int start_vert)
{
    MSTParentArray parents(n);

    struct mstprim_node
    {
        float weight;
        int vert;

        bool operator<(const mstprim_node& other) const
        {
            return other.weight < weight;
        }
    };

    std::priority_queue<mstprim_node> heap;
    std::vector<float> weights(n, std::numeric_limits<float>::max());
    FlagArr in_mst(n);

    // Initialize 
    heap.push({0.0f, start_vert});
    weights[start_vert] = 0.0f;

    // Insert a new vertex into the MST at each iteration,
    // when all vertices are inserted, do n iterations
    for (int i = 0; i < n; i++)
    {
        // Pick the vertex with minimum weight
        const auto [v_weight, v] = heap.top();
        heap.pop();
        in_mst.Set(v);

        // Update neighbours (== every vertex in the graph)
        for (int u = 0; u < n; u++)
        {
            if (u == v || in_mst.Query(u))
                continue;

            float& u_weight = weights[u];
            float new_uweight = distance_matrix.get(v, u);

            if (u_weight > new_uweight)
            {
                u_weight = new_uweight;
                parents[u] = v;
                heap.push({u_weight, u});
            }
        }
    }

    return parents;
}

MSTEdgeArr BuildMSTEdgeArr(int n, const MSTParentArray& parents)
{
    MSTEdgeArr edge_arr;
    edge_arr.reserve(n);

    for (int i = 0; i < n; i++)
    {
        const int parent = parents[i];

        if (i == parent)
            continue;

        edge_arr.emplace_back(MST_Edge{i, parent});
        edge_arr.emplace_back(MST_Edge{parent, i});
    }

    std::sort(edge_arr.begin(), edge_arr.end());
    return edge_arr;
}

int GetFirstAdjVert(const MSTEdgeArr& edges, int v)
{
    const auto it = std::lower_bound(edges.begin(), edges.end(), MST_Edge{v, 0});
    assert(it != edges.end());
    return it - edges.begin();
}