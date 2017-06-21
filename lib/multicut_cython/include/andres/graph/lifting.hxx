#pragma once
#ifndef ANDRES_GRAPH_LIFTING_HXX
#define ANDRES_GRAPH_LIFTING_HXX

#include <cassert>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <iterator> // std::iterator_traits
#include <algorithm> // std::fill
#include <vector>

#include "grid-graph.hxx"
#include "bfs.hxx"


namespace andres {
namespace graph {

enum class LiftingMetric { PathLength, L2 };

/// Lift a graph.
template<class INPUT_GRAPH, class OUTPUT_GRAPH>
inline void
lift(
    const INPUT_GRAPH& inputGraph,
    OUTPUT_GRAPH& outputGraph,
    const std::size_t distanceUpperBound,
    const std::size_t distanceLowerBound = 0
) {
    typedef std::size_t size_type;

    if(outputGraph.numberOfVertices() != 0)
        throw std::runtime_error("output graph is not empty.");
    
    outputGraph.insertVertices(inputGraph.numberOfVertices());

    BreadthFirstSearchData<size_type> breadthFirstSearchData(inputGraph.numberOfVertices());
    std::vector<std::size_t> visited;
    std::vector<std::size_t> vertices;

    for (size_type v = 0; v < inputGraph.numberOfVertices(); ++v)
    {
        breadthFirstSearch(
            inputGraph,
            v,
            [&](size_type w, size_type depth, bool& proceed, bool& add)
            {
                proceed = true;
                add = false;
                
                if (depth <= distanceUpperBound)
                {
                    if (depth + 1 <= distanceUpperBound)
                    {
                        add = true;
                        visited.push_back(w);        
                    }

                    if (depth > distanceLowerBound)
                        vertices.push_back(w);
                }
            },
            breadthFirstSearchData
            );

        std::sort(vertices.begin(), vertices.end());

        for (auto w : vertices)
            outputGraph.insertEdge(v, w);

        vertices.clear();

        breadthFirstSearchData.depth(v) = BreadthFirstSearchData<size_type>::NOT_VISITED;
        for (auto w : visited)
            breadthFirstSearchData.depth(w) = BreadthFirstSearchData<size_type>::NOT_VISITED;
        
        visited.clear();
    }
}

/// Lift a grid graph - a faster implementation using the grid structure.
template<class INPUT_GRAPH_VISITOR, class OUTPUT_GRAPH>
inline void
lift(
    const GridGraph<2, INPUT_GRAPH_VISITOR>& inputGraph,
    OUTPUT_GRAPH& outputGraph,
    const std::size_t distanceUpperBound,
    const std::size_t distanceLowerBound = 0,
    LiftingMetric metric = LiftingMetric::PathLength
) {
    typedef GridGraph<2, INPUT_GRAPH_VISITOR> INPUT_GRAPH;

    typedef std::size_t size_type;
    typedef typename INPUT_GRAPH::VertexCoordinate VertexCoordinate;

    const size_type distanceUpperBoundSquared = distanceUpperBound * distanceUpperBound;
    const size_type distanceLowerBoundSquared = distanceLowerBound * distanceLowerBound;

    if(outputGraph.numberOfVertices() != 0)
        throw std::runtime_error("output graph is not empty.");

    outputGraph.insertVertices(inputGraph.numberOfVertices());

    VertexCoordinate cv;
    for (size_type v = 0; v < inputGraph.numberOfVertices(); ++v)
    {
        inputGraph.vertex(v, cv);

        // fill above portion of the window
        if (cv[1] > 0)
        {
            const std::size_t row0 = cv[1] < distanceUpperBound ? 0 : cv[1] - distanceUpperBound;
            std::size_t offsetY = 1;

            // We use yPlus = y+1 to avoid the 0-1 case. In this block y = yPlus-1.
            for(std::size_t yPlus = cv[1]; yPlus > row0; --yPlus, ++offsetY)
            {
                const std::size_t offsetX = (metric == LiftingMetric::PathLength) ?
                                            distanceUpperBound - offsetY:
                                            ::floor(::sqrt(distanceUpperBoundSquared - offsetY * offsetY));

                const std::size_t col0 = cv[0] < offsetX ? 0 : cv[0] - offsetX;
                std::size_t colN = cv[0] + offsetX;
                if(colN > inputGraph.shape(0) - 1)
                    colN = inputGraph.shape(0) - 1;
                
                for (std::size_t x = col0; x <= colN; ++x)
                {
                    if (metric == LiftingMetric::PathLength)
                    {
                        const std::size_t distance = ::abs(x - cv[0]) + ::abs(yPlus - 1 - cv[1]);

                        if (distance > distanceLowerBound)
                        {
                            const size_type w = inputGraph.vertex({{x, yPlus - 1}});
                            outputGraph.insertEdge(v, w);
                        }
                    }
                    else
                    {
                        const std::size_t sqaredDistance = (x - cv[0]) * (x - cv[0]) + (yPlus - 1 - cv[1]) * (yPlus - 1 - cv[1]);

                        if (sqaredDistance > distanceLowerBoundSquared)
                        {
                            const size_type w = inputGraph.vertex({{x, yPlus - 1}});
                            outputGraph.insertEdge(v, w);
                        }
                    }
                }
            }
        }

        // Add middle horizontal line, except for the central point
        {
            const std::size_t offsetX = distanceUpperBound;
            const std::size_t y = cv[1];
            const std::size_t col0 = cv[0] < offsetX ? 0 : cv[0] - offsetX;
            std::size_t colN = cv[0] + offsetX;

            if (colN > inputGraph.shape(0) - 1)
                colN = inputGraph.shape(0) - 1;
            
            if (cv[0] > distanceLowerBound)
                for (std::size_t x = col0; x <= cv[0] - distanceLowerBound - 1; ++x)
                {
                    const size_type& w = inputGraph.vertex({{x, y}});
                    outputGraph.insertEdge(v, w);
                }

            for (std::size_t x = cv[0] + distanceLowerBound + 1; x <= colN; ++x)
            {
                const size_type& w = inputGraph.vertex({{x, y}});
                outputGraph.insertEdge(v, w);
            }
        }

        // fill below window
        if (cv[1] < inputGraph.shape(1) - 1)
        {
            const std::size_t row0 = cv[1] + 1;
            std::size_t rowN = cv[1] + distanceUpperBound;
            
            if (cv[1] + distanceUpperBound > inputGraph.shape(1) - 1)
                rowN = inputGraph.shape(1) - 1;

            std::size_t offsetY = 1;
            for (std::size_t y = row0; y <= rowN; ++y, ++offsetY)
            {
                const std::size_t offsetX = (metric == LiftingMetric::PathLength) ?
                                            distanceUpperBound - offsetY :
                                            ::floor(::sqrt(distanceUpperBoundSquared - offsetY * offsetY));
                
                const std::size_t col0 = cv[0] < offsetX ? 0 : cv[0] - offsetX;
                std::size_t colN = cv[0] + offsetX;

                if (colN > inputGraph.shape(0) - 1)
                    colN = inputGraph.shape(0) - 1;

                for (std::size_t x = col0; x <= colN; ++x)
                {
                    if (metric == LiftingMetric::PathLength)
                    {
                        const std::size_t distance = ::abs(x - cv[0]) + ::abs(y - cv[1]);
                        
                        if (distance > distanceLowerBound)
                        {
                            const size_type w = inputGraph.vertex({{x, y}});
                            outputGraph.insertEdge(v, w);
                        }
                    }
                    else
                    {
                        const std::size_t sqaredDistance = (x - cv[0]) * (x - cv[0]) + (y - cv[1]) * (y - cv[1]);

                        if (sqaredDistance > distanceLowerBoundSquared)
                        {
                            const size_type w = inputGraph.vertex({{x, y}});
                            outputGraph.insertEdge(v, w);
                        }
                    }
                }
            }
        }
    }
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_LIFTING_HXX
