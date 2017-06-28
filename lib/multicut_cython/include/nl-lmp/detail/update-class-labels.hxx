#pragma once
#ifndef NL_LMP_UPDATE_CLASS_LABELS_HXX
#define NL_LMP_UPDATE_CLASS_LABELS_HXX

#include <queue>
#include <vector>

#include "../problem.hxx"
#include "../solution.hxx"



namespace nl_lmp
{

namespace detail
{

template<typename GRAPH>
inline
double update_class_labels(Problem<GRAPH> const& problem, Solution const& input, Solution& output)
{
    struct entry
    {
        entry(double __gain, size_t __v, size_t __k, size_t __edition = 0) :
            gain(__gain), v(__v), k(__k), edition(__edition)
        {}
        
        bool operator <(entry const& other) const
        {
            return gain < other.gain;
        }

        double gain;
        size_t v;
        size_t k;
        size_t edition;
    };

    auto const number_of_labels = problem.numberOfClasses();

    output = input;

    std::priority_queue<entry> Q;
    std::vector<double> gains(problem.numberOfVertices() * number_of_labels);

    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
    {
        double best_gain = .0;
        size_t best_label = 0;

        auto const v_label = input[v].classIndex;

        for (size_t k = 0; k < number_of_labels; ++k)
        {
            auto gain = problem.getUnaryCost(v, v_label) - problem.getUnaryCost(v, k);

            for (auto it = problem.liftedGraph().adjacenciesFromVertexBegin(v); it != problem.liftedGraph().adjacenciesFromVertexEnd(v); ++it)
            {
                auto const e = it->edge();
                auto const w = it->vertex();

                auto const w_label = input[w].classIndex;

                if (input[v].clusterIndex == input[w].clusterIndex)
                    gain += problem.getPairwiseJoinCost(v, w, v_label, w_label) - problem.getPairwiseJoinCost(v, w, k, w_label);
                else
                    gain += problem.getPairwiseCutCost(v, w, v_label, w_label) - problem.getPairwiseCutCost(v, w, k, w_label);
            }

            if (gain > best_gain)
            {
                best_gain = gain;
                best_label = k;
            }

            gains[v*number_of_labels + k] = gain;
        }

        if (best_gain > 1e-6)
            Q.push(entry(best_gain, v, best_label));
    }

    std::vector<char> visited(problem.numberOfVertices());
    std::vector<size_t> vertex_editions(problem.numberOfVertices());

    double total_gain = .0;
    while (!Q.empty())
    {
        auto const ent = Q.top();
        Q.pop();

        if (visited[ent.v] || ent.edition < vertex_editions[ent.v])
            continue;

        total_gain += ent.gain;

        auto const v_former_label = output[ent.v].classIndex;
        output[ent.v].classIndex = ent.k;
        visited[ent.v] = 1;
        
        for (auto it = problem.liftedGraph().adjacenciesFromVertexBegin(ent.v); it != problem.liftedGraph().adjacenciesFromVertexEnd(ent.v); ++it)
        {
            auto const w = it->vertex();

            if (!visited[w])
            {
                auto const e = it->edge();
                auto const w_label = output[w].classIndex;

                double best_gain = .0;
                size_t best_label = 0;

                for (size_t k = 0; k < number_of_labels; ++k)
                {
                    // cancel previous contribution and add the new one
                    if (output[ent.v].clusterIndex == output[w].clusterIndex)
                    {
                        gains[w*number_of_labels + k] -= problem.getPairwiseJoinCost(w, ent.v, w_label, v_former_label, e) - problem.getPairwiseJoinCost(w, ent.v, k, v_former_label, e);
                        gains[w*number_of_labels + k] += problem.getPairwiseJoinCost(w, ent.v, w_label, ent.k, e) - problem.getPairwiseJoinCost(w, ent.v, k, ent.k, e);
                    }
                    else
                    {
                        gains[w*number_of_labels + k] -= problem.getPairwiseCutCost(w, ent.v, w_label, v_former_label, e) - problem.getPairwiseCutCost(w, ent.v, k, v_former_label, e);
                        gains[w*number_of_labels + k] += problem.getPairwiseCutCost(w, ent.v, w_label, ent.k, e) - problem.getPairwiseCutCost(w, ent.v, k, ent.k, e);
                    }

                    if (gains[w*number_of_labels + k] > best_gain)
                    {
                        best_gain = gains[w*number_of_labels + k];
                        best_label = k;
                    }
                }

                ++vertex_editions[w];

                if (best_gain > 1e-6)
                    Q.push(entry(best_gain, w, best_label, vertex_editions[w]));
            }
        }
    }

    return total_gain;
}

} // of namespace detail

} // of namespace nl_lmp

#endif
