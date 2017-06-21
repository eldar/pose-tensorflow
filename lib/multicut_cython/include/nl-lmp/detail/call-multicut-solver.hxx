#pragma once
#ifndef NL_LMP_CALL_MULTICUT_SOLVER_HXX
#define NL_LMP_CALL_MULTICUT_SOLVER_HXX

#include <andres/graph/complete-graph.hxx>
#include <andres/graph/multicut-lifted/kernighan-lin.hxx>
#include <andres/graph/multicut/kernighan-lin.hxx>

#include "compute-objective.hxx"



namespace nl_lmp
{

namespace detail
{

template<typename GRAPH, typename ECA, typename VERTEXLABELS>
inline
std::vector<size_t> call_kernighanLin(Problem<GRAPH> const& problem, ECA const& edge_weights, VERTEXLABELS const& vertex_cluster_labels)
{
    return andres::graph::multicut_lifted::kernighanLin(problem.originalGraph(), problem.liftedGraph(), edge_weights, vertex_cluster_labels);
}

template<typename GRAPHVISITOR, typename ECA, typename VERTEXLABELS>
inline
std::vector<size_t> call_kernighanLin(Problem<andres::graph::CompleteGraph<GRAPHVISITOR>> const& problem, ECA const& edge_weights, VERTEXLABELS const& vertex_cluster_labels)
{
    return andres::graph::multicut::kernighanLin(problem.originalGraph(), edge_weights, vertex_cluster_labels);
}

template<typename GRAPH>
inline
double update_multicut(Problem<GRAPH> const& problem, Solution const& input, Solution& output)
{
    // given class labels frozen we're essentially optimizing the multicut objective
    // therefore, I just make use of existing Kernighan-Lin code
    // this is also much faster, because accessing 'edge_weights' is way faster

    auto energy = compute_obj_value(problem, input);

    std::vector<size_t> vertex_cluster_labels(problem.numberOfVertices());
    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
    {
        vertex_cluster_labels[v] = input[v].clusterIndex;

        output[v].classIndex = input[v].classIndex;
    }

    std::vector<double> edge_weights(problem.liftedGraph().numberOfEdges());
    for (size_t e = 0; e < problem.liftedGraph().numberOfEdges(); ++e)
    {
        auto const v0 = problem.liftedGraph().vertexOfEdge(e, 0);
        auto const v1 = problem.liftedGraph().vertexOfEdge(e, 1);

        edge_weights[e] = problem.getPairwiseCutCost(v0, v1, input[v0].classIndex, input[v1].classIndex, e) - problem.getPairwiseJoinCost(v0, v1, input[v0].classIndex, input[v1].classIndex, e);
    }

    auto new_vertex_labels = call_kernighanLin(problem, edge_weights, vertex_cluster_labels);

    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
        output[v].clusterIndex = new_vertex_labels[v];

    double energy_after = compute_obj_value(problem, output);

    return energy - energy_after;
}

}

}

#endif
