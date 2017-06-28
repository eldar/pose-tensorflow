#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_ILP_CALLBACK_HXX
#define ANDRES_GRAPH_MULTICUT_ILP_CALLBACK_HXX

#include <vector>
#include <deque>
#include <array>
#include <algorithm> // std::copy

#include "andres/graph/complete-graph.hxx"
#include "andres/graph/components.hxx"
#include "andres/graph/paths.hxx"
#include "andres/graph/shortest-paths.hxx"


namespace andres {
namespace graph {
namespace multicut
{

/// Solver for the Minimum Cost Multicut Problem for arbitrary graphs.
///
/// This is a variant of the solver proposed in 
/// 
/// Andres B., Kroeger T., Briggman K. L., Denk W., Korogod N., Knott G., Koethe U. and Hamprecht F. A.
/// Globally Optimal Closed-surface Segmentation for Connectomics. ECCV 2012
/// http://dx.doi.org/10.1007/978-3-642-33712-3_56
///
/// This code operates on graphs whereas the code used to produce the results 
/// in the above publication operates on cellular complexes. While cellular 
/// complexes are richer structures than graphs and facilitates additional 
/// tweaks of the solver (e.g. cell suppression), graphs are more common. 
/// This code has a wider range of applications.
/// 
template<typename ILP, typename GRAPH, typename ECA, typename ELA>
inline
void ilp_callback(GRAPH const& graph, ECA const& edgeCosts, ELA const& inputLabels, ELA& outputLabels, size_t timeLimitSeconds = 86400)
{
    struct Visitor
    {
        bool operator()(ELA const& edge_labels) const
        {
            return true;
        }
    } visitor;

    ilp_callback<ILP>(graph, edgeCosts, inputLabels, outputLabels, visitor, timeLimitSeconds);
}

template<typename ILP, typename GRAPH, typename ECA, typename ELA, typename VIS>
inline
void ilp_callback(GRAPH const& graph, ECA const& edgeCosts, ELA const& inputLabels, ELA& outputLabels, VIS& visitor, size_t timeLimitSeconds = 86400)
{
    struct SubgraphWithCut
    {
        SubgraphWithCut(ILP const& ilp)
            : ilp_(ilp) 
            {}

        bool vertex(size_t v) const
        {
            return true;
        }

        bool edge(size_t e) const
        {
            return ilp_.label(e) < .5;
        }

        ILP const& ilp_;
    };

    class Callback: public ILP::Callback
    {
    public:
        Callback(ILP& solver, GRAPH const& graph) :
            ILP::Callback(solver), graph_(graph)
        {}

        void separateAndAddLazyConstraints() override
        {
            ComponentsBySearch<GRAPH> components;
            std::deque<size_t> path;
            std::vector<ptrdiff_t> buffer;
            std::vector<double> variables(graph_.numberOfEdges());
            std::vector<double> coefficients(graph_.numberOfEdges());

            components.build(graph_, SubgraphWithCut(*this));

            // search for violated non-chordal cycles and add corresp. inequalities
            for (size_t edge = 0; edge < graph_.numberOfEdges(); ++edge) 
                if (this->label(edge) > .5)
                {
                    auto v0 = graph_.vertexOfEdge(edge, 0);
                    auto v1 = graph_.vertexOfEdge(edge, 1);

                    if (components.areConnected(v0, v1))
                    { 
                        // search for shortest path
                        spsp(graph_, SubgraphWithCut(*this), v0, v1, path, buffer);
                        
                        // skip chordal paths
                        if (findChord(graph_, path.begin(), path.end(), true).first)
                            continue;

                        // add inequality
                        auto sz = path.size();

                        for (size_t j = 0; j < sz - 1; ++j)
                        {
                            variables[j] = static_cast<double>(graph_.findEdge(path[j], path[j + 1]).second);
                            coefficients[j] = 1.0;
                        }

                        variables[sz-1] = static_cast<double>(edge);
                        coefficients[sz-1] = -1.0;

                        this->addLazyConstraint(variables.begin(), variables.begin() + sz, coefficients.begin(), 0, std::numeric_limits<double>::infinity());
                    }
                }
        }
    private:
        struct SubgraphWithCut
        {
            SubgraphWithCut(Callback& callback)
                : callback_(callback) 
                {}

            bool vertex(size_t v) const
            {
                return true;
            }

            bool edge(size_t e) const
            {
                return callback_.label(e) < .5;
            }

            Callback& callback_;
        };

        GRAPH const& graph_;
    };

    ILP ilp;

    ilp.setRelativeGap(0.0);
    ilp.setAbsoluteGap(0.0);
    ilp.setTimeLimit(timeLimitSeconds);
    ilp.addVariables(edgeCosts.size(), edgeCosts.data());

    Callback callback(ilp, graph);
    ilp.setCallback(callback);

    ilp.optimize();

    ComponentsBySearch<GRAPH> components;
    components.build(graph, SubgraphWithCut(ilp));

    for (size_t edge = 0; edge < graph.numberOfEdges(); ++edge)
    {
        auto v0 = graph.vertexOfEdge(edge, 0);
        auto v1 = graph.vertexOfEdge(edge, 1);

        outputLabels[edge] = components.areConnected(v0, v1) ? 0 : 1;
    }
}



/// Solver for the Minimum Cost Multicut Problem for complete graphs (Set Partition Problem).
template<typename ILP, typename GRAPH_VISITOR, typename ECA, typename ELA>
inline
void ilp_callback(CompleteGraph<GRAPH_VISITOR> const& graph, ECA const& edgeCosts, ELA const& inputLabels, ELA& outputLabels, size_t timeLimitSeconds = 86400)
{
    struct Visitor
    {
        bool operator()(ELA const& edge_labels) const
        {
            return true;
        }
    } visitor;

    ilp_callback<ILP>(graph, edgeCosts, inputLabels, outputLabels, visitor);
}

template<typename ILP, typename GRAPH_VISITOR, typename ECA, typename ELA, typename VIS>
inline
void ilp_callback(CompleteGraph<GRAPH_VISITOR> const& graph, ECA const& edgeCosts, ELA const& inputLabels, ELA& outputLabels, VIS& visitor, size_t timeLimitSeconds = 86400)
{
    struct SubgraphWithCut
    {
        SubgraphWithCut(ILP const& ilp)
            : ilp_(ilp) 
            {}

        bool vertex(size_t v) const
        {
            return true;
        }

        bool edge(size_t e) const
        {
            return ilp_.label(e) < .5;
        }

        ILP const& ilp_;
    };

    class Callback: public ILP::Callback
    {
    public:
        Callback(ILP& solver, CompleteGraph<GRAPH_VISITOR> const& graph) :
            ILP::Callback(solver), graph_(graph)
        {}

        void separateAndAddLazyConstraints() override
        {
            std::array<double, 3> variables;
            std::array<double, 3> coefficients;

            for (size_t edge = 0; edge < graph_.numberOfEdges(); ++edge) 
                if (this->label(edge) > .5)
                {
                    variables[2] = edge;

                    auto v0 = graph_.vertexOfEdge(edge, 0);
                    auto v1 = graph_.vertexOfEdge(edge, 1);

                    for (size_t i = 0; i < graph_.numberOfVertices(); ++i)
                    {
                        if (i == v0 || i == v1)
                            continue;

                        variables[0] = graph_.findEdge(v0, i).second;
                        variables[1] = graph_.findEdge(v1, i).second;

                        if (this->label(variables[0]) < .5 && this->label(variables[1]) < .5)
                        {
                            coefficients[0] =  1.0;
                            coefficients[1] =  1.0;
                            coefficients[2] = -1.0;

                            this->addLazyConstraint(variables.begin(), variables.end(), coefficients.begin(), 0, std::numeric_limits<double>::infinity());
                        }
                    }
                }
        }
    private:
        CompleteGraph<GRAPH_VISITOR> const& graph_;
    };

    ILP ilp;

    ilp.setRelativeGap(0.0);
    ilp.setAbsoluteGap(0.0);
    ilp.setTimeLimit(timeLimitSeconds);
    ilp.addVariables(edgeCosts.size(), edgeCosts.data());

    Callback callback(ilp, graph);
    ilp.setCallback(callback);

    ilp.optimize();

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

#endif // #ifndef ANDRES_GRAPH_MULTICUT_ILP_CALLBACK_HXX
