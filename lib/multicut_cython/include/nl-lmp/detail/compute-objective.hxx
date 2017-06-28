#pragma once
#ifndef NL_LMP_COMPUTE_OBJECTIVE_HXX
#define NL_LMP_COMPUTE_OBJECTIVE_HXX

#include "../problem.hxx"
#include "../solution.hxx"



namespace nl_lmp
{

namespace detail
{

template<typename GRAPH>
inline
double compute_obj_value(Problem<GRAPH> const& problem, Solution const& solution)
{
    double obj_value = .0;

    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
        obj_value += problem.getUnaryCost(v, solution[v].classIndex);

    for (size_t e = 0; e < problem.liftedGraph().numberOfEdges(); ++e)
    {
        auto const v0 = problem.liftedGraph().vertexOfEdge(e, 0);
        auto const v1 = problem.liftedGraph().vertexOfEdge(e, 1);

        if (solution[v0].clusterIndex != solution[v1].clusterIndex)
            obj_value += problem.getPairwiseCutCost(v0, v1, solution[v0].classIndex, solution[v1].classIndex);
        else
            obj_value += problem.getPairwiseJoinCost(v0, v1, solution[v0].classIndex, solution[v1].classIndex);
    }

    return obj_value;
}

}

}

#endif 
