#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_ILP_HXX
#define ANDRES_GRAPH_MULTICUT_ILP_HXX

#include <vector>
#include <deque>
#include <array>
#include <algorithm>

#include <andres/graph/complete-graph.hxx>
#include <andres/graph/components.hxx>
#include <andres/graph/paths.hxx>
#include <andres/graph/shortest-paths.hxx>

namespace andres {
namespace graph {
namespace multicut {

/// Algorithm for the Minimum Cost Multicut Problem.
///
/// A variant of the algorithm proposed in
///
/// Andres B., Kroeger T., Briggman K. L., Denk W., Korogod N., Knott G.,
/// Koethe U. and Hamprecht F. A. Globally Optimal Closed-surface Segmentation
/// for Connectomics. ECCV 2012. http://dx.doi.org/10.1007/978-3-642-33712-3_56
///
template<typename ILP, typename GRAPH, typename ECA, typename ELA>
inline void
ilp(
    GRAPH const& graph,
    ECA const& edgeCosts,
    ELA const& inputLabels,
    ELA& outputLabels,
    size_t numberOfIterations = std::numeric_limits<size_t>::max()
) {
    struct Visitor {
        bool operator()(ELA const& edge_labels) const {
            return true;
        }
    } visitor;

    ilp<ILP>(graph, edgeCosts, inputLabels, outputLabels, visitor, numberOfIterations);
}

template<typename ILP, typename GRAPH, typename ECA, typename ELA, typename VIS>
inline void
ilp(
    GRAPH const& graph,
    ECA const& edgeCosts,
    ELA const& inputLabels,
    ELA& outputLabels,
    VIS& visitor,
    size_t numberOfIterations = std::numeric_limits<size_t>::max()
) {
    struct SubgraphWithCut {
        SubgraphWithCut(const ILP& ilp)
            : ilp_(ilp) 
            {}
        bool vertex(const size_t v) const
            { return true; }
        bool edge(const size_t e) const
            { return ilp_.label(e) < .5; }
        const ILP& ilp_;
    };

    ComponentsBySearch<GRAPH> components;
    ILP ilp;
    std::deque<size_t> path;
    std::vector<ptrdiff_t> buffer;
    std::vector<size_t> variables(graph.numberOfEdges());
    std::vector<double> coefficients(graph.numberOfEdges());

    auto addCycleInequalities = [&] ()
    {
        components.build(graph, SubgraphWithCut(ilp));

        // search for violated non-chordal cycles and add corresp. inequalities
        size_t nCycle = 0;

        for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge) 
            if (ilp.label(edge) > .5)
            {
                auto v0 = graph.vertexOfEdge(edge, 0);
                auto v1 = graph.vertexOfEdge(edge, 1);

                if (components.areConnected(v0, v1))
                { 
                    // search for shortest path
                    spsp(graph, SubgraphWithCut(ilp), v0, v1, path, buffer);
                    
                    // skip chordal paths
                    if (findChord(graph, path.begin(), path.end(), true).first)
                        continue;

                    // add inequality
                    auto sz = path.size();

                    for (size_t j = 0; j < path.size() - 1; ++j)
                    {
                        variables[j] = graph.findEdge(path[j], path[j + 1]).second;
                        coefficients[j] = 1.0;
                    }

                    variables[path.size() - 1] = edge;
                    coefficients[path.size() - 1] = -1.0;

                    ilp.addConstraint(variables.begin(), variables.begin() + path.size(), coefficients.begin(), 0, std::numeric_limits<double>::infinity());

                    ++nCycle;
                }
            }

        return nCycle;
    };

    auto repairSolution = [&] ()
    {
        for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge)
        {
            auto v0 = graph.vertexOfEdge(edge, 0);
            auto v1 = graph.vertexOfEdge(edge, 1);

            outputLabels[edge] = components.areConnected(v0, v1) ? 0 : 1;
        }

        ilp.setStart(outputLabels.begin());
    };

    ilp.initModel(graph.numberOfEdges(), edgeCosts.data());
    ilp.setStart(inputLabels.begin());

    for (size_t i = 0; numberOfIterations == 0 || i < numberOfIterations; ++i)
    {
        if (i != 0)
        {
            repairSolution();

            if (!visitor(outputLabels))
                break;
        }

        ilp.optimize();

        if (addCycleInequalities() == 0)
            break;
    }

    repairSolution();
}

/// Algorithm for the Set Partition Problem.
///
/// The Set Partition Problem is the Minimum Cost Multicut Problem for complete
/// graphs.
///
template<typename ILP, typename GRAPH_VISITOR, typename ECA, typename ELA>
inline void
ilp(
    CompleteGraph<GRAPH_VISITOR> const& graph,
    ECA const& edgeCosts,
    ELA const& inputLabels,
    ELA& outputLabels,
    size_t numberOfIterations = std::numeric_limits<size_t>::max()
) {
    struct Visitor
    {
        bool operator()() const
        {
            return true;
        }
    } visitor;

    ilp<ILP>(graph, edgeCosts, inputLabels, outputLabels, visitor, numberOfIterations);
}

template<typename ILP, typename GRAPH_VISITOR, typename ECA, typename ELA, typename VIS>
inline void
ilp(
    CompleteGraph<GRAPH_VISITOR> const& graph,
    ECA const& edgeCosts,
    ELA const& inputLabels,
    ELA& outputLabels,
    VIS& visitor,
    size_t numberOfIterations = std::numeric_limits<size_t>::max()
) {
    struct SubgraphWithCut {
        SubgraphWithCut(const ILP& ilp)
            : ilp_(ilp) 
            {}
        bool vertex(const size_t v) const
            { return true; }
        bool edge(const size_t e) const
            { return ilp_.label(e) < .5; }
        const ILP& ilp_;
    };

    ILP ilp;
    std::array<double, 3> variables;
    std::array<double, 3> coefficients;

    auto addCycleInequalities = [&] ()
    {
        size_t nCycle = 0;

        for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge) 
            if (ilp.label(edge) > .5)
            {
                variables[2] = edge;

                auto v0 = graph.vertexOfEdge(edge, 0);
                auto v1 = graph.vertexOfEdge(edge, 1);

                for (size_t i = 0; i < graph.numberOfVertices(); ++i)
                {
                    if (i == v0 || i == v1)
                        continue;

                    variables[0] = graph.findEdge(v0, i).second;
                    variables[1] = graph.findEdge(v1, i).second;

                    if (ilp.label(variables[0]) < .5 && ilp.label(variables[1]) < .5)
                    {
                        coefficients[0] =  1.0;
                        coefficients[1] =  1.0;
                        coefficients[2] = -1.0;

                        ilp.addConstraint(variables.begin(), variables.end(), coefficients.begin(), 0, std::numeric_limits<double>::infinity());

                        ++nCycle;
                    }
                }
            }

        return nCycle;
    };

    ilp.initModel(graph.numberOfEdges(), edgeCosts.data());
    ilp.setStart(inputLabels.begin());

    for (size_t i = 0; numberOfIterations == 0 || i < numberOfIterations; ++i)
    {
        if (!visitor())
            break;

        ilp.optimize();

        if (addCycleInequalities() == 0)
            break;
    }

    ComponentsBySearch<CompleteGraph<GRAPH_VISITOR>> components;
    components.build(graph, SubgraphWithCut(ilp));

    for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge)
    {
        auto v0 = graph.vertexOfEdge(edge, 0);
        auto v1 = graph.vertexOfEdge(edge, 1);

        outputLabels[edge] = components.areConnected(v0, v1) ? 0 : 1;
    }
}

} // namespace multicut
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_ILP_HXX
