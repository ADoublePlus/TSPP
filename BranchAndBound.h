#pragma once

#include "Flag_Array.h"
#include "MST.h"
#include "Matrix.h"

#include <algorithm>
#include <queue>
#include <stack>
#include <vector>

namespace branch_and_bound
{
    /*
        Search tree starts with an empty set, and at each step, each node branches into two nodes, one that includes 'e' and one that doesn't.
        'e' is selected in the increasing length order.
        Node is discarded if the lower bound for the given node is greater than the current best solution or approximation found.
     */
    enum Constraint_Types : uint8_t
    {
        // Let the default initialized value be Allowed
        Allowed,
        Disallowed,
        Required,
        _NUM_CONSTRAINTS
    };

    // Reserve 2 bits for each constraint 
    static_assert(_NUM_CONSTRAINTS < 4);

    typedef int edge_t;

    inline edge_t EdgeId(int n, int row, int col)
    {
        return row * n + col;
    }

    struct Vertex_Pair
    {
        int row;
        int col;
    };

    inline Vertex_Pair EdgeId(int n, edge_t edge_id)
    {
        std::div_t res = std::div(edge_id, n);
        return {res.quot, res.rem};
    }

    struct Edge
    {
        edge_t edge_id;
        float weight;

        bool operator<(const Edge& other) const
        {
            return weight > other.weight;
        }
    };

    enum TSPP_Type
    {
        FixedStartVert,
        ArbitraryStartVert
    };

    // Holds reserved data for the constraints for all edges (n^2)
    typedef Bitset<2> Constraints;

    // Instead of carrying n^2 number of constraints, carry the diff from the parent
    // Tightening is applied later on
    struct Node
    {
        int8_t row;
        int8_t col;
        int8_t diff_constraint;
        int8_t counter;
        Constraints constraints;
    };

    static_assert((1 << 8) > MAX_CITIES);

    /*
        Tighten constraints are given these three rules:
            - When all but two edges adjacent to a vertex are excluded, those two edges have to be included as otherwise it would be possible for a tour to exist.
            - When two edges adjacent to a vertex are included, all other edges adjacent to this vertex have to be excluded.
            - When including an edge would complete a subtour with other included edges, this edge has to be excluded.
        
        Note that a Constraint_Diff takes an 'Allowed' edge and marks it whatever constraint member is.
        Reversing a patch is as simple as marking edge as 'Allowed' again.
     */
    struct Constraint_Diff
    {
        int8_t row;
        int8_t col;
        int16_t constraint;

        bool operator<(const Constraint_Diff& other) const
        {
            const auto [_row, _col] = std::minmax(row, col);
            const auto [_oRow, _oCol] = std::minmax(other.row, other.col);

            if (_row == _oRow)
                return _col < _oCol;

            return _row < _oRow;
        }

        bool operator==(const Constraint_Diff& other) const
        {
            return (row == other.row && col == other.col) || (row == other.col && col == other.row);
        }
    };

    void PrintConstraints(const int n, const Constraints& constraints);
    inline void ApplyConstraint(const int n, Constraints* constraints, const Constraint_Diff& diff);
    inline void RevertConstraint(const int n, Constraints* constraints, const Constraint_Diff& diff);

    typedef std::vector<Constraint_Diff> ConstraintTighteningPatch;
    void ApplyConstraintTighteningPatch(const int n, Constraints* constraints,
                                        const ConstraintTighteningPatch& patch);
    void RevertConstraintTighteningPatch(const int n, Constraints* constraints,
                                         const ConstraintTighteningPatch& patch);

    enum TightenConstraints_ReturnCode
    {
        Valid,
        ValidAndLeaf,
        NonValid
    };

    // Cycle checks are done during the traversal phase
    template <TSPP_Type tspp_type>
    TightenConstraints_ReturnCode IterativelyTightenConstraints(const int n, Constraints* constraints,
                                                                ConstraintTighteningPatch* patch,
                                                                const int start_vert);

    template <TSPP_Type tspp_type>
    TightenConstraints_ReturnCode TightenConstraints(const int n, Constraints* constraints,
                                                     ConstraintTighteningPatch* patch,
                                                     const int start_vert);

    template TightenConstraints_ReturnCode TightenConstraints<TSPP_Type::ArbitraryStartVert>(const int n, const Constraints& constraints,
                                                                                             ConstraintTighteningPatch* patch,
                                                                                             const int start_vert);

    template TightenConstraints_ReturnCode TightenConstraints<TSPP_Type::FixedStartVert>(const int n, const Constraints& constraints,
                                                                                         ConstraintTighteningPatch* patch,
                                                                                         const int start_vert);

    typedef std::vector<Constraint_Types> DBGExpConstraint;
    DBGExpConstraint DBExpandConstraints(const int n, const Constraints& constraints);

    float CalculateSolutionLength(const int n, const Constraints& constraints, const Matrix& distance_matrix);

    // Find the shortest length edge e that is allowed 
    Vertex_Pair CandidateEdge(const int n, const Constraints& constraints, const Matrix& distance_matrix);

    template <TSPP_Type tspp_type>
    float LowerBound(const Constraints& node, const Matrix& distance_matrix, const int start_vert);
    template float LowerBound<TSPP_Type::ArbitraryStartVert>(const Constraints& node, const Matrix& distance_matrix, const int start_vert);
    template float LowerBound<TSPP_Type::FixedStartVert>(const Constraints& node, const Matrix& distance_matrix, const int start_vert);

    // Solve the TSPP problem for a complete graph with given distance matrix
    template <TSPP_Type tspp_type>
    Path SolveMetricTSPP(const int n, const Matrix& distance_matrix,
                         const int start_vert, const int exec_time_limit);

    template Path SolveMetricTSPP<TSPP_Type::ArbitraryStartVert>(const int n, const Matrix& distance_matrix,
                                                                           const int start_vert, const int exec_time_limit);

    template Path SolveMetricTSPP<TSPP_Type::FixedStartVert>(const int n, const Matrix& distance_matrix,
                                                             const int start_vert, const int exec_time_limit);

    // Bidirectional search, if a path is found over required edges, insertion of the new edge (u, v) will create a sub-tour 
    bool YieldsCycle(const int n, const Constraints& constraints, int u, int v);

    template <TSPP_Type tspp_type>
    void AssertConstraintsValidity(const int n, const Constraints& constraints, const int start_vert);
    template void AssertConstraintsValidity<TSPP_Type::ArbitraryStartVert>(const int n, const Constraints& constraints, const int start_vert);
    template void AssertConstraintsValidity<TSPP_Type::FixedStartVert>(const int n, const Constraints& constraints, const int start_vert);

    Path ConstraintsToPath(const int n, const Constraints& constraints, const int start_vert);
} // Namespace branch_and_bound