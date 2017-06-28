#pragma once
#ifndef ANDRES_GRAPH_GRAPH_HDF5_HXX
#define ANDRES_GRAPH_GRAPH_HDF5_HXX

#include <stdexcept>
#include <string>
#include <vector>

#include "hdf5.hxx"
#include "../graph.hxx"

namespace andres {
namespace graph {
namespace hdf5 {

template<>
template<class VISITOR>
struct GraphTraitsHDF5<Graph<VISITOR> > {
    static const int ID;
};
template<class VISITOR>
    const int GraphTraitsHDF5<Graph<VISITOR> >::ID = 10000;

template <class VISITOR>
void save(const hid_t, const std::string&, const Graph<VISITOR>&);

template <class VISITOR>
void load(const hid_t, const std::string&, Graph<VISITOR>&);

template <class VISITOR>
void
save(
    const hid_t parentHandle,
    const std::string& graphName,
    const Graph<VISITOR>& graph
) {
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(parentHandle, graphName, true);

    try {
        save(groupHandle, "graph-type-id", GraphTraitsHDF5<Graph<VISITOR> >::ID);
        save(groupHandle, "multiple-edges-enabled", static_cast<unsigned char>(graph.multipleEdgesEnabled()));
        save(groupHandle, "number-of-vertices", graph.numberOfVertices());
        save(groupHandle, "number-of-edges", graph.numberOfEdges());
        if(graph.numberOfEdges() != 0) {
            std::vector<std::size_t> vecIJ;
            vecIJ.resize(2 * graph.numberOfEdges());
            std::size_t *ptrI = &vecIJ[0];
            std::size_t *ptrJ = &vecIJ[graph.numberOfEdges()];
            for(std::size_t e=0; e<graph.numberOfEdges(); ++e) {
                *(ptrI++) = graph.vertexOfEdge(e, 0);
                *(ptrJ++) = graph.vertexOfEdge(e, 1);
            }
            save(groupHandle, "edges", {graph.numberOfEdges(), 2}, &vecIJ[0]);
        }
    } catch (std::exception& e) {
        closeGroup(groupHandle);
        throw std::runtime_error("error saving graph: " + std::string(e.what()));
    }

    closeGroup(groupHandle);
}

template <class VISITOR>
void
load(
    const hid_t parentHandle,
    const std::string& graphName,
    Graph<VISITOR>& graph
) {
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(parentHandle, graphName);

    std::string sError;

    try {
        int id = 0;
        load(groupHandle, "graph-type-id", id);
        if(id != GraphTraitsHDF5<Graph<VISITOR> >::ID) {
            sError = "graph type id mismatch.";
            goto cleanup;
        }

        std::size_t numberOfVertices = 0;
        load(groupHandle, "number-of-vertices", numberOfVertices);

        std::size_t numberOfEdges = 0;
        load(groupHandle, "number-of-edges", numberOfEdges);

        graph.assign(numberOfVertices);

        if(numberOfEdges != 0) {
            graph.reserveEdges(numberOfEdges);

            std::vector<std::size_t> bufferIJ;
            {
                std::vector<std::size_t> shape;
                load(groupHandle, "edges", shape, bufferIJ);
                if(shape.size() != 2) {
                    sError = "edges dataset is not 2-dimensional.";
                    goto cleanup;
                }
                if(shape[0] != numberOfEdges || shape[1] != 2) {
                    sError = "edges dataset has incorrect shape.";
                    goto cleanup;
                }
                assert(shape[0] * 2 == bufferIJ.size());
            }
            std::size_t* ptrI = &bufferIJ[0];
            std::size_t* ptrJ = &bufferIJ[numberOfEdges];
            {
                unsigned char multipleEdgesEnabled;
                load(groupHandle, "multiple-edges-enabled", multipleEdgesEnabled);
                graph.multipleEdgesEnabled() = static_cast<bool>(multipleEdgesEnabled);
            }
            for(std::size_t i=0;i<numberOfEdges;++i) {
                const std::size_t s = *(ptrI++);
                const std::size_t t = *(ptrJ++);
                if(s >= numberOfVertices || t >= numberOfVertices) {
                    sError = "vertex index out of bounds.";
                    goto cleanup;
                }
                const std::size_t e = graph.insertEdge(s, t);
            }
        }
    } catch(std::exception& e) {
        sError = e.what();
    }

cleanup:
    closeGroup(groupHandle);
    if(!sError.empty()) {
        throw std::runtime_error("error loading graph: " + sError);
    }
}

} //namespace hdf
} //namespace graph
} //namespace andres


#endif // #ifndef ANDRES_GRAPH_GRAPH_HDF5_HXX
