#pragma once
#ifndef NL_LMP_ILP_CALLBACK_HXX
#define NL_LMP_ILP_CALLBACK_HXX

#include <andres/graph/components.hxx>
#include <andres/graph/paths.hxx>
#include <andres/graph/shortest-paths.hxx>

#include "problem.hxx"
#include "solution.hxx"



namespace nl_lmp
{

template<typename ILP, typename GRAPH, typename VISITOR>
inline
Solution ilp_callback(Problem<GRAPH> const& problem, Solution const& input, size_t time_limit_seconds = 43200)
{
    struct EmptyVisitor
    {
    } visitor;

    return ilp_callback<ILP>(problem, input, visitor, time_limit_seconds);
}

template<typename ILP, typename GRAPH, typename VISITOR>
inline
Solution ilp_callback(Problem<GRAPH> const& problem, Solution const& input, VISITOR& visitor, size_t time_limit_seconds = 43200)
{
    struct SubgraphWithCut
    {
        SubgraphWithCut(ILP const& ilp, size_t __offsetY) :
            ilp_(ilp), offsetY(__offsetY)
        {}

        bool vertex(size_t v) const
        {
            return true;
        }

        bool edge(size_t e) const
        {
            return ilp_.label(offsetY + e) < .5;
        }

        ILP const& ilp_;
        size_t const offsetY;
    };

    size_t const offsetY = problem.numberOfVertices()*problem.numberOfClasses();
    size_t const offsetQ = offsetY + problem.originalGraph().numberOfEdges();
    size_t const offsetZ = offsetQ + problem.originalGraph().numberOfEdges()*problem.numberOfClasses()*problem.numberOfClasses();

    auto x = [&problem](size_t v, size_t l)
    {
        return v*problem.numberOfClasses() + l;
    };
    
    auto y = [&problem, &offsetY](size_t v0, size_t v1)
    {
        return offsetY + problem.originalGraph().findEdge(v0, v1).second;
    };
    
    auto q = [&problem, &offsetQ](size_t v0, size_t v1, size_t c0, size_t c1)
    {
        if (v1 < v0)
            std::swap(c0, c1);

        return offsetQ + problem.originalGraph().findEdge(v0, v1).second*problem.numberOfClasses()*problem.numberOfClasses() + c0*problem.numberOfClasses() + c1;
    };

    auto z = [&problem, &offsetZ](size_t v0, size_t v1, size_t c0, size_t c1)
    {
        if (v1 < v0)
            std::swap(c0, c1);

        return offsetZ + problem.originalGraph().findEdge(v0, v1).second*problem.numberOfClasses()*problem.numberOfClasses() + c0*problem.numberOfClasses() + c1;
    };

    class Callback: public ILP::Callback
    {
    public:
        Callback(ILP& solver, Problem<GRAPH> const& problem, VISITOR& visitor) :
            ILP::Callback(solver), problem_(problem),
            offsetY_(problem.numberOfVertices()*problem.numberOfClasses()),
            visitor_(visitor)
        {

        }

        void separateAndAddLazyConstraints() override
        {
            std::cerr << '\t' << this->objectiveBound_
                << '\t' << this->objectiveBest_
                << '\t' << this->objective()
                << std::flush;

            size_t n = 0;

            andres::graph::ComponentsBySearch<GRAPH> components;
            components.build(problem_.originalGraph(), SubgraphWithCut(*this, offsetY_));

            std::deque<size_t> path;
            std::vector<ptrdiff_t> buffer;
            std::vector<size_t> variables(problem_.originalGraph().numberOfEdges());
            std::vector<double> coefficients(problem_.originalGraph().numberOfEdges());

            for (size_t edge = 0; edge < problem_.originalGraph().numberOfEdges(); ++edge)
            {
                auto const v0 = problem_.originalGraph().vertexOfEdge(edge, 0);
                auto const v1 = problem_.originalGraph().vertexOfEdge(edge, 1);

                if (this->label(y(v0, v1)) > .5 && components.areConnected(v0, v1))
                {
                    // search for shortest path
                    spsp(problem_.originalGraph(), SubgraphWithCut(*this, offsetY_), v0, v1, path, buffer);
                    
                    // skip chordal paths
                    if (findChord(problem_.originalGraph(), path.begin(), path.end(), true).first)
                        continue;

                    // add inequality
                    for (size_t j = 0; j < path.size() - 1; ++j)
                    {
                        variables[j] = y(path[j], path[j + 1]);
                        coefficients[j] = 1.0;
                    }

                    variables[path.size() - 1] = y(v0, v1);
                    coefficients[path.size() - 1] = -1.0;

                    this->addLazyConstraint(variables.begin(), variables.begin() + path.size(), coefficients.begin(), 0, std::numeric_limits<double>::infinity());

                    ++n;
                }
            }

            // the following is optimization for the complete graph
            // uncomment this one and comment the previous for-loop

            // std::array<double, 3> const coefficients({ 1.0, 1.0, -1.0 });

            // std::array<size_t, 3> variables;
            // for (size_t edge = 0; edge < problem_.originalGraph().numberOfEdges(); ++edge)
            // {
            //     auto const v0 = problem_.originalGraph().vertexOfEdge(edge, 0);
            //     auto const v1 = problem_.originalGraph().vertexOfEdge(edge, 1);

            //     if (this->label(y(v0, v1)) > .5)
            //     {
            //         variables[2] = y(v0, v1);

            //         for (size_t i = 0; i < problem_.originalGraph().numberOfVertices(); ++i)
            //         {
            //             if (i == v0 || i == v1)
            //                 continue;

            //             variables[0] = y(i, v0);
            //             variables[1] = y(i, v1);

            //             if (this->label(y(i, v0)) < .5 && this->label(y(i, v1)) < .5)
            //             {
            //                 this->addLazyConstraint(variables.begin(), variables.end(), coefficients.begin(), .0, std::numeric_limits<double>::infinity());
            //                 ++n;
            //             }
            //         }
            //     }
            // }

            std::cerr << '\t' << n << std::endl;

            if (n == 0)
            {
                andres::graph::ComponentsBySearch<GRAPH> components;
                components.build(problem_.originalGraph(), SubgraphWithCut(*this, offsetY_));

                Solution solution(problem_.numberOfVertices());

                for (size_t v = 0; v < problem_.numberOfVertices(); ++v)
                {
                    solution[v].clusterIndex = components.labels_[v];
                    solution[v].classIndex = problem_.numberOfClasses();

                    for (size_t c = 0; c < problem_.numberOfClasses(); ++c)
                        if (this->label(v*problem_.numberOfClasses() + c) > .5)
                        {
                            solution[v].classIndex = c;                
                            break;
                        }
                }
            }
        }
    private:
        struct SubgraphWithCut
        {
            SubgraphWithCut(Callback& callback, size_t __offsetY) :
                callback_(callback), offsetY(__offsetY)
            {}

            bool vertex(size_t v) const
            {
                return true;
            }

            bool edge(size_t e) const
            {
                return callback_.label(offsetY + e) < .5;
            }

            Callback& callback_;
            size_t const offsetY;
        };

        size_t y(size_t v0, size_t v1) const
        {
            return offsetY_ + problem_.originalGraph().findEdge(v0, v1).second;
        }

        Problem<GRAPH> const& problem_;
        size_t const offsetY_;

        VISITOR& visitor_;
    };


    ILP ilp;

    ilp.setRelativeGap(0.0);
    ilp.setAbsoluteGap(0.0);
    ilp.setTimeLimit(time_limit_seconds);

    {
        std::vector<double> coefficients(offsetZ * 2 - offsetQ);

        // set coefficients of variables x
        for (size_t v = 0; v < problem.numberOfVertices(); ++v)
            for (size_t l = 0; l < problem.numberOfClasses(); ++l)
                coefficients[x(v, l)] = problem.getUnaryCost(v, l);

        // there are no coefficients for variables y, e.g. they are zero
        // the values of y are derived from z using constraints

        // set coefficients of variables q and z
        for (size_t edge = 0; edge < problem.originalGraph().numberOfEdges(); ++edge)
        {
            auto const v0 = problem.originalGraph().vertexOfEdge(edge, 0);
            auto const v1 = problem.originalGraph().vertexOfEdge(edge, 1);
        
            for (size_t c0 = 0; c0 < problem.numberOfClasses(); ++c0)
                for (size_t c1 = 0; c1 < problem.numberOfClasses(); ++c1)
                {
                    coefficients[q(v0, v1, c0, c1)] = problem.getPairwiseJoinCost(v0, v1, c0, c1);

                    coefficients[z(v0, v1, c0, c1)] = problem.getPairwiseCutCost(v0, v1, c0, c1) - problem.getPairwiseJoinCost(v0, v1, c0, c1);
                }
        }

        ilp.addVariables(coefficients.size(), coefficients.data());

        std::cout   << "   " << offsetY                         << " variables x" << std::endl
                    << "   " << offsetQ - offsetY               << " variables y" << std::endl
                    << "   " << offsetZ - offsetQ               << " variables q" << std::endl
                    << "   " << coefficients.size() - offsetZ   << " variables z" << std::endl;
    }

    {
        std::cout << "adding all impossible part class constraints: " << std::flush;

        size_t n = 0;

        size_t variables[1];
        double const  coefficients[] = { 1.0 };
        for (size_t v = 0; v < problem.numberOfVertices(); ++v)
            for (size_t l = 0; l < problem.numberOfClasses(); ++l)
                if (problem.getUnaryCost(v, l) > 500.0)
                {
                    variables[0] = x(v, l);
                    
                    ilp.addConstraint(variables, variables + 1, coefficients, .0, .0);
                    ++n;
                }

        std::cout << n << std::endl;
    }

    bool mustSelect = true;

    if (!mustSelect)
    {
        std::cout << "adding all coupling constraints: " << std::flush;

        size_t n = 0;

        std::vector<size_t> variables(1 + problem.numberOfClasses());
        std::vector<double> coefficients(1 + problem.numberOfClasses(), 1.0);

        coefficients.back() = -1.0;

        for (size_t edge = 0; edge < problem.originalGraph().numberOfEdges(); ++edge)
        {
            auto const v0 = problem.originalGraph().vertexOfEdge(edge, 0);
            auto const v1 = problem.originalGraph().vertexOfEdge(edge, 1);

            variables.back() = y(v0, v1);

            for (size_t c = 0; c < problem.numberOfClasses(); ++c)
                variables[c] = x(v0, c);

            ilp.addConstraint(variables.begin(), variables.end(), coefficients.begin(), .0, std::numeric_limits<double>::infinity());
            ++n;
            
            for (size_t c = 0; c < problem.numberOfClasses(); ++c)
                variables[c] = x(v1, c);

            ilp.addConstraint(variables.begin(), variables.end(), coefficients.begin(), .0, std::numeric_limits<double>::infinity());
            ++n;
        }

        std::cout << n << std::endl;
    }

    {
        std::cout << "adding all uniqueness constraints: " << std::flush;

        size_t n = 0;

        vector<size_t> variables(problem.numberOfClasses());
        vector<size_t> const coefficients(problem.numberOfClasses(), 1.0);

        for (size_t v = 0; v < problem.numberOfVertices(); ++v)
        {
            for (size_t l = 0; l < problem.numberOfClasses(); ++l)
                variables[l] = x(v, l);

            ilp.addConstraint(variables.begin(), variables.end(), coefficients.begin(), mustSelect ? 1.0 : .0, 1.0);
            ++n;
        }

        std::cout << n << std::endl;
    }

    {
        std::cout << "adding all linearization constraints: " << std::flush;

        size_t n = 0;

        size_t vi[] = { 0, 0, 0, 0 };

        for (size_t edge = 0; edge < problem.originalGraph().numberOfEdges(); ++edge)
        {
            auto const v0 = problem.originalGraph().vertexOfEdge(edge, 0);
            auto const v1 = problem.originalGraph().vertexOfEdge(edge, 1);

            vi[0] = y(v0, v1);

            {
                double const coefficients[] = { 1.0, 1.0, 1.0, -1.0 };

                for (size_t c0 = 0; c0 < problem.numberOfClasses(); ++c0)
                {
                    vi[1] = x(v0, c0);

                    for (size_t c1 = 0; c1 < problem.numberOfClasses(); ++c1)
                    {
                        vi[2] = x(v1, c1);
                        vi[3] = z(v0, v1, c0, c1);

                        ilp.addConstraint(vi, vi + 4, coefficients, -std::numeric_limits<double>::infinity(), 2.0);
                        ++n;

                        vi[3] = q(v0, v1, c0, c1);

                        ilp.addConstraint(vi + 1, vi + 4, coefficients + 1, -std::numeric_limits<double>::infinity(), 1.0);
                        ++n;
                    }
                }
            }

            {
                double const coefficients[] = {-1.0, 1.0};

                for (size_t c0 = 0; c0 < problem.numberOfClasses(); ++c0)
                    for (size_t c1 = 0; c1 < problem.numberOfClasses(); ++c1)
                    {
                        vi[1] = z(v0, v1, c0, c1);

                        ilp.addConstraint(vi, vi + 2, coefficients, -std::numeric_limits<double>::infinity(), .0);
                        ++n;
                    }
            }

            {
                double const coefficients[] = { 1.0, -1.0 };

                for (size_t c0 = 0; c0 < problem.numberOfClasses(); ++c0)
                    for (size_t c1 = 0; c1 < problem.numberOfClasses(); ++c1)
                    {
                        vi[0] = q(v0, v1, c0, c1);

                        vi[1] = x(v0, c0);

                        ilp.addConstraint(vi, vi + 2, coefficients, -std::numeric_limits<double>::infinity(), .0);
                        ++n;

                        vi[1] = x(v1, c1);

                        ilp.addConstraint(vi, vi + 2, coefficients, -std::numeric_limits<double>::infinity(), .0);
                        ++n;


                        vi[0] = z(v0, v1, c0, c1);

                        vi[1] = x(v0, c0);

                        ilp.addConstraint(vi, vi + 2, coefficients, -std::numeric_limits<double>::infinity(), .0);
                        ++n;

                        vi[1] = x(v1, c1);

                        ilp.addConstraint(vi, vi + 2, coefficients, -std::numeric_limits<double>::infinity(), .0);
                        ++n;
                    }
            }
        }

        std::cout << n << std::endl;
    }

    Callback callback(ilp, problem, visitor);
    ilp.setCallback(callback);

    ilp.optimize();

    andres::graph::ComponentsBySearch<GRAPH> components;
    components.build(problem.originalGraph(), SubgraphWithCut(ilp, offsetY));

    Solution output(problem.numberOfVertices());

    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
    {
        output[v].clusterIndex = components.labels_[v];
        output[v].classIndex = problem.numberOfClasses();

        for (size_t c = 0; c < problem.numberOfClasses(); ++c)
            if (ilp.label(x(v, c)) > .5)
            {
                output[v].classIndex = c;
                break;
            }
    }

    return output;
}

}

#endif
