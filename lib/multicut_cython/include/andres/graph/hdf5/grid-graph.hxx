#pragma once
#ifndef ANDRES_GRAPH_GRID_GRAPH_HDF5_HXX
#define ANDRES_GRAPH_GRID_GRAPH_HDF5_HXX

#include <stdexcept>
#include <string>
#include <array>
#include <vector>

#include "hdf5.hxx"
#include "../grid-graph.hxx"

namespace andres {
namespace graph {
namespace hdf5 {

template<>
template<unsigned char D, class VISITOR>
struct GraphTraitsHDF5<GridGraph<D, VISITOR> > {
    static const int ID;
};
template<unsigned char D, class VISITOR>
    const int GraphTraitsHDF5<GridGraph<D, VISITOR> >::ID = 10003;

template<unsigned char D, class VISITOR>
void save(const hid_t, const std::string&, const GridGraph<D, VISITOR>&);

template<unsigned char D, class VISITOR>
void load(const hid_t, const std::string&, GridGraph<D, VISITOR>&);

template<unsigned char D, class VISITOR>
void
save(
    const hid_t parentHandle,
    const std::string& graphName,
    const GridGraph<D, VISITOR>& graph
) {
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(parentHandle, graphName, true);

    std::string sError;
    try {
        save(groupHandle, "graph-type-id", GraphTraitsHDF5<GridGraph<D, VISITOR> >::ID);

        std::array<std::size_t, D> shape;
        for(std::size_t i=0; i<D; ++i) {
            shape[i] = graph.shape(i);
        }
        save(groupHandle, "shape", {D}, &shape[0]);
    } catch(std::exception& e) {
        sError = e.what();
    }
    closeGroup(groupHandle);
    if(!sError.empty()) {
        throw std::runtime_error("error saving grid graph: " + sError);
    }
}

template<unsigned char D, class VISITOR>
void
load(
    const hid_t parentHandle,
    const std::string& graphName,
    GridGraph<D, VISITOR>& graph
) {
    HandleCheck<ANDRES_GRAPH_HDF5_DEBUG> handleCheck;
    hid_t groupHandle = openGroup(parentHandle, graphName);
    
    std::string sError;
    try {
        int id;
        load(groupHandle, "graph-type-id", id);
        if(id != GraphTraitsHDF5<GridGraph<D, VISITOR> >::ID) {
            sError = "graph type id mismatch.";
            goto cleanup;
        }

        std::vector<std::size_t> nDims;
        std::vector<std::size_t> shape;
        load(groupHandle, "shape", nDims, shape);
        if(nDims.size() != 1 || nDims[0] != D) {
            sError = "shape mismatch.";
            goto cleanup;
        }

        typename GridGraph<D, VISITOR>::VertexCoordinate vc;
        std::copy(shape.begin(), shape.end(), vc.begin());
        graph.assign(vc);
    } catch(std::exception& e) {
        sError = "failed to load grid graph: " + std::string(e.what());
        goto cleanup;
    }

cleanup:
    closeGroup(groupHandle);
    if(!sError.empty()) {
        throw std::runtime_error("error loading grid graph: " + sError);
    }
}

} //namespace hdf
} //namespace graph
} //namespace andres

#endif // #ifndef ANDRES_GRAPH_GRID_GRAPH_HDF5_HXX
