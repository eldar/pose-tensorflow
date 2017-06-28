#pragma once

#include <andres/marray.hxx>

// ILP is currently solving only pure multicut (not lifted and without class labels)
enum class SolvingMethod
{ 
    Simple,
    Joint
#ifdef WITH_GUROBI
    , ILP
#endif
};

// if unaries.shape(1) == 1, then it only clusters
andres::Marray<size_t> solve_nl_lmp_complete_graph(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, SolvingMethod method, bool do_logit_transform, bool with_suppression = false);

andres::Marray<size_t> solve_nl_lmp_sparse_graph(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, SolvingMethod method, bool do_logit_transform, bool with_suppression = false);

// with partial solution
andres::Marray<size_t> solve_nl_lmp_complete_graph(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, andres::View<size_t> const& partial_solution, SolvingMethod method, bool do_logit_transform, bool with_suppression = false);

andres::Marray<size_t> solve_nl_lmp_sparse_graph(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, andres::View<size_t> const& partial_solution, SolvingMethod method, bool do_logit_transform, bool with_suppression = false);


double compute_objective(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, andres::View<size_t> const& solution, bool do_logit_transform);