#include "Approx.h"
#include "Flag_Array.h"
#include "MST.h"

#include <iostream>
#include <stack>
#include <vector>

namespace mst_approx
{
    using namespace data;

    Path Approximate(const int num_verts, const Matrix& distance_matrix, const int start_vert)
    {
        MSTParentArray mst_parents = MSTPrim(num_verts, distance_matrix, start_vert);
        MSTEdgeArr edge_arr = BuildMSTEdgeArr(num_verts, mst_parents);
        const int num_edges = 2 * num_verts - 2;
        assert(edge_arr.size() == num_edges);

        FlagArr is_visited(num_verts);
        std::stack<int> stack;
        stack.push(start_vert);
        Path ts_path;
        ts_path.reserve(num_verts);

        while (!stack.empty())
        {
            const int u = stack.top();
            stack.pop();
            is_visited.Set(u);
            ts_path.push_back(u);

            int adj = GetFirstAdjVert(edge_arr, u);

            for (int e = adj; e < num_edges; e++)
            {
                MST_Edge& edge = edge_arr[e];

                if (edge.u != u)
                    break;

                int v = edge.v;

                if (!is_visited.Query(v))
                {
                    stack.push(v);
                }
            }
        }

        return ts_path;
    }
} // Namespace mst_approx

namespace iterative_improv
{
    void improv_2opt(const int num_verts, const Matrix& distance_matrix,
                     const Opt_Strategy strategy, Path* path)
    {
        const int num_edges = num_verts - 1;

        bool improv_found = true;
        float total_gain = 0.0f;

        printf("2-OPT w/ MST 2-approx as input\n");
        float initial_length = geometry::CalculatePathLength(*path, distance_matrix);
        printf("Initial length: %f\n", initial_length);

        float max_gain = 0.0f;

        struct Opt2_Strategy
        {
            int i;
            int j;
        } max_gain_segment;

        const float epsilon = 1e-4;

        // Iterate until no local improvement is found 
        int step = 0;

        while (improv_found)
        {
            improv_found = false;
            max_gain = epsilon;

            for (int i = 0; i < num_edges; i++)
            {
                const int e1u = (*path)[i];
                const int e1v = (*path)[i + 1];
                float edge1_len = distance_matrix.get(e1u, e1v);

                for (int j = i + 2; j < num_edges; j++)
                {
                    const int e2u = (*path)[j];
                    const int e2v = (*path)[j + 1];
                    float edge2_len = distance_matrix.get(e2u, e2v);

                    float gain = edge1_len + edge2_len - distance_matrix.get(e1u, e2u) - 
                                                         distance_matrix.get(e1v, e2v);

                    // Do the first move that has positive gain
                    if (strategy == Opt_Strategy::FirstGain)
                    {
                        if (gain > epsilon)
                        {
                            improv_found = true;
                            float _l0 = geometry::CalculatePathLength(*path, distance_matrix);
                            std::reverse(path->begin() + i + 1, path->begin() + j + 1);
                            float _l1 = geometry::CalculatePathLength(*path, distance_matrix);
                            assert(fabs(_l0 - gain - _l1) < epsilon);
                            total_gain += gain;

                            break;
                        }
                    }
                    else if (strategy == Opt_Strategy::BestGain)
                    {
                        if (gain > max_gain)
                        {
                            max_gain = gain;
                            max_gain_segment = {i, j};
                        }
                    }
                }
                
                if (strategy == Opt_Strategy::FirstGain && improv_found)
                    break;
            }

            if (strategy == Opt_Strategy::BestGain)
            {
                // Apply the best opt2 move
                std::reverse(path->begin() + max_gain_segment.i + 1,
                             path->begin() + max_gain_segment.j + 1);

                total_gain += max_gain;
            }

            step++;
        }

        float _l = geometry::CalculatePathLength(*path, distance_matrix);
        printf("2OPT, Final length: %f | %f\n", initial_length - total_gain, _l);
    }

    void improv_3opt(const int num_verts, const Matrix& distance_matrix,
                     const Opt_Strategy strategy, Path* path)
    {
        const int num_edges = num_verts - 1;

        bool improv_found = true;
        float total_gain = 0.0f;

        printf("3-OPT w/ MST 2-approx as input\n");
        float initial_length = geometry::CalculatePathLength(*path, distance_matrix);
        printf("Initial length: %f\n", initial_length);

        float max_gain = 0.0f;

        enum Opt3_Move_Case
        {
            case_0, // No change
            case_1, // One possible 2-opt move, one segment left unchanged
            case_2, // One possible 2-opt move, one segment left unchanged
            case_3, // One possible 2-opt move, one segment left unchanged
            case_4 // Equivalent to three or more 2-opt moves. ad, be, cf
        };

        struct Opt3_Segment
        {
            int i;
            int j;
            int k;

            Opt3_Move_Case move_case;
        } max_gain_segment;

        const float epsilon = 1e-4;

        auto ApplyOpt3Move = [](Path* _path, Opt3_Segment opt3_move)
        {
            if (opt3_move.move_case == Opt3_Move_Case::case_0) {}
            else if (opt3_move.move_case == Opt3_Move_Case::case_1)
            {
                std::reverse(_path->begin() + opt3_move.i + 1, _path->begin() + opt3_move.j + 1);
            }
            else if (opt3_move.move_case == Opt3_Move_Case::case_2)
            {
                std::reverse(_path->begin() + opt3_move.j + 1, _path->begin() + opt3_move.k + 1);
            }
            else if (opt3_move.move_case == Opt3_Move_Case::case_3)
            {
                std::reverse(_path->begin() + opt3_move.i + 1, _path->begin() + opt3_move.k + 1);
            }
            else if (opt3_move.move_case == Opt3_Move_Case::case_4)
            {
                Path _cpy(opt3_move.k - opt3_move.i);

                for (int x = 0, _end = opt3_move.k - opt3_move.j; x < _end; x++)
                {
                    _cpy[x] = (*_path)[opt3_move.j + 1 + x];
                }

                for (int x = 0, _end = opt3_move.j - opt3_move.i; x < _end; x++)
                {
                    _cpy[opt3_move.k - opt3_move.j + x] = (*_path)[opt3_move.i + 1 + x];
                }

                for (int x = 0, _end = opt3_move.k - opt3_move.i; x < _end; x++)
                {
                    (*_path)[opt3_move.i + 1 + x] = _cpy[x];
                }
            }
        };

        // Iterate until no local improvement is found
        int step = 0;

        while (improv_found)
        {
            improv_found = false;
            max_gain = epsilon;

            for (int i = 0; i < num_edges; i++)
            {
                const int A = (*path)[i];
                const int B = (*path)[i + 1];
                float d_AB = distance_matrix.get(A, B);

                for (int j = i + 2; j < num_edges; j++)
                {
                    const int C = (*path)[j];
                    const int D = (*path)[j + 1];
                    float d_CD = distance_matrix.get(C, D);
                    float d_AC = distance_matrix.get(A, C);
                    float d_AD = distance_matrix.get(A, D);
                    float d_BD = distance_matrix.get(B, D);

                    for (int k = j + 2; k < num_edges; k++)
                    {
                        const int E = (*path)[k];
                        const int F = (*path)[k + 1];
                        float d_EF = distance_matrix.get(E, F);
                        float d_EB = distance_matrix.get(E, B);
                        float d_FB = distance_matrix.get(F, B);
                        float d_CE = distance_matrix.get(C, E);
                        float d_DF = distance_matrix.get(D, F);
                        float d_EA = distance_matrix.get(E, A);
                        float d_CF = distance_matrix.get(C, F);

                        float d0 = d_AB + d_CD + d_EF;
                        float d1 = d_AC + d_BD + d_EF;
                        float d2 = d_AB + d_CE + d_DF;
                        float d3 = d_AD + d_EB + d_CF;
                        float d4 = d_FB + d_CD + d_EA;

                        Opt3_Move_Case c = Opt3_Move_Case::case_0;
                        float gain = 0.0f;

                        if (d0 > d1)
                        {
                            c = Opt3_Move_Case::case_1;
                            gain = d0 - d1;
                        }
                        else if (d0 > d2)
                        {
                            c = Opt3_Move_Case::case_2;
                            gain = d0 - d2;
                        }
                        else if (d0 > d4)
                        {
                            c = Opt3_Move_Case::case_3;
                            gain = d0 - d4;
                        }
                        else if (d0 > d3)
                        {
                            c = Opt3_Move_Case::case_4;
                            gain = d0 - d3;
                        }

                        // Do first move that has positive gain
                        if (strategy == Opt_Strategy::FirstGain)
                        {
                            if (gain > epsilon)
                            {
                                improv_found = true;
                                float _l0 = geometry::CalculatePathLength(*path, distance_matrix);
                                Path _path_cpy(path->begin(), path->end());
                                ApplyOpt3Move(path, {i, j, k, c});
                                float _l1 = geometry::CalculatePathLength(*path, distance_matrix);
                                assert(fabs(_l0 - gain - _l1) < epsilon);
                                total_gain += gain;

                                break;
                            }
                        }
                        else if (strategy == Opt_Strategy::BestGain)
                        {
                            if (gain > max_gain)
                            {
                                max_gain = gain;
                                max_gain_segment = {i, k, j, c};
                            }
                        }
                    }

                    if (strategy == Opt_Strategy::FirstGain && improv_found)
                        break;
                }

                if (strategy == Opt_Strategy::FirstGain && improv_found)
                    break;
            }

            if (strategy == Opt_Strategy::BestGain)
            {
                // Apply the best opt3 move
                ApplyOpt3Move(path, max_gain_segment);
                total_gain += max_gain;
            }

            step++;
        }

        float _l = geometry::CalculatePathLength(*path, distance_matrix);
        printf("3OPT, Final length: %f | %f\n", initial_length - total_gain, _l);
    }
} // Namespace iterative_improv