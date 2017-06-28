#pragma once
#ifndef NL_LMP_GREEDY_ADDITIVE_HXX
#define NL_LMP_GREEDY_ADDITIVE_HXX

#include <vector>

#include <andres/graph/multicut-lifted/greedy-additive.hxx>

#include "problem.hxx"
#include "solution.hxx"



namespace nl_lmp
{

// this function computes only multicut, so it expects already initialized class labeling
template<typename GRAPH>
inline
Solution greedyAdditiveEdgeContraction(Problem<GRAPH> const& problem, Solution const& input_labeling)
{
    std::vector<double> edge_weights(problem.liftedGraph().numberOfEdges());
    for (size_t e = 0; e < problem.liftedGraph().numberOfEdges(); ++e)
    {
        auto const v0 = problem.liftedGraph().vertexOfEdge(e, 0);
        auto const v1 = problem.liftedGraph().vertexOfEdge(e, 1);

        edge_weights[e] = problem.getPairwiseCutCost(v0, v1, input_labeling[v0].classIndex, input_labeling[v1].classIndex, e) - problem.getPairwiseJoinCost(v0, v1, input_labeling[v0].classIndex, input_labeling[v1].classIndex, e);
    }

    auto vertex_cluster_labels = andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(problem.originalGraph(), problem.liftedGraph(), edge_weights);

    Solution solution(problem.numberOfVertices());
    for (size_t i = 0; i < problem.numberOfVertices(); ++i)
        solution[i].clusterIndex = vertex_cluster_labels[i];

    return solution;
}

}

#endif 
