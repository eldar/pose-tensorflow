#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LP_HXX
#define ANDRES_GRAPH_MULTICUT_LP_HXX

#include <vector>
#include <deque>
#include <algorithm>

#include <andres/graph/complete-graph.hxx>
#include <andres/graph/shortest-paths.hxx>


namespace andres {
namespace graph {
namespace multicut {

template<typename LP, typename GRAPH, typename ECA>
inline
std::vector<double> lp(GRAPH const& graph, ECA const& edgeCosts, size_t numberOfIterations = std::numeric_limits<size_t>::max())
{
    struct Visitor
    {
        bool operator()() const
        {
            return true;
        }
    } visitor;

    return lp<LP>(graph, edgeCosts, visitor, numberOfIterations);
}

template<typename LP, typename GRAPH, typename ECA, typename VIS>
inline
std::vector<double> lp(GRAPH const& graph, ECA const& edgeCosts, VIS& visitor, size_t numberOfIterations = std::numeric_limits<size_t>::max())
{
    const double tolerance = std::numeric_limits<float>::epsilon();

    LP lp;

    std::vector<double> coefficients(graph.numberOfEdges());
    std::vector<double> distances(graph.numberOfVertices());
    std::vector<size_t> parents(graph.numberOfVertices());
    std::deque<size_t> path;
    std::vector<size_t> variables(graph.numberOfEdges());
    std::vector<double> vars(graph.numberOfEdges());

    auto addCycleInequalities = [&] ()
    {
        for (size_t i = 0; i < vars.size(); ++i)
            vars[i] = std::min(std::max(.0, lp.variableValue(i)), 1.0);
            // although in Gurobi we constrain the variables to be in the [0,1] range,
            // sometimes Gurobi finds sligthly negative solutions of the order of 1e-13.
            // The latter totally screws up Dijkstra's shortest path algorithm

        size_t nCycle = 0;
        for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge) 
        {
            auto v0 = graph.vertexOfEdge(edge, 0);
            auto v1 = graph.vertexOfEdge(edge, 1);

            // search for shortest path
            double distance;
            spsp(graph, DefaultSubgraphMask<>(), v0, v1, vars.begin(), path, distance, distances.begin(), parents.begin());

            bool chordless = true;
            for (auto it1 = path.begin(); it1 != path.end() - 2 && chordless; ++it1)
                for (auto it2 = it1 + 2; it2 != path.end(); ++it2)
                {
                    if (it1 == path.begin() && it2 == path.end() - 1)
                        continue;

                    auto e = graph.findEdge(*it1, *it2);
                    if (e.first && std::min(std::max(.0, lp.variableValue(e.second)), 1.0) > distances[*it2] - distances[*it1] + tolerance)
                    {
                        chordless = false;
                        break;
                    }
                }

            if (!chordless)
                continue;

            if (vars[edge] > distance + tolerance)
            {
                // add inequality
                for (size_t j = 0; j < path.size() - 1; ++j)
                {
                    coefficients[j] = 1.0;
                    variables[j] = graph.findEdge(path[j], path[j + 1]).second;
                }

                coefficients[path.size() - 1] = -1.0;
                variables[path.size() - 1] = edge;

                lp.addConstraint(variables.begin(), variables.begin() + path.size(), coefficients.begin(), .0, std::numeric_limits<double>::infinity());

                ++nCycle;
            }
        }

        return nCycle;
    };

    lp.initModel(graph.numberOfEdges(), edgeCosts.data());

    for (size_t i = 0; numberOfIterations == 0 || i < numberOfIterations; ++i)
    {
        if (!visitor())
            break;

        lp.optimize();

        if (addCycleInequalities() == 0)
            break;
    }

    std::vector<double> edge_values(graph.numberOfEdges());

    for (size_t i = 0; i < graph.numberOfEdges(); ++i)
        edge_values[i] = lp.variableValue(i);

    return edge_values;
}

}
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_LP_HXX
