#include <stdexcept>

#include <andres/functional.hxx>
#include <andres/graph/complete-graph.hxx>
#include <andres/graph/dfs.hxx>
#include <andres/graph/graph.hxx>
#ifdef WITH_GUROBI
    #include <andres/ilp/gurobi-callback.hxx>
    #include <andres/graph/multicut/ilp-callback.hxx>
#endif


#include <nl-lmp/greedy-additive.hxx>
#include <nl-lmp/solve-alternating.hxx>
#include <nl-lmp/solve-joint.hxx>

#include "nl-lmp.hxx"


using namespace std;

template<typename GRAPH, typename VLA, typename ELA>
inline
void edgeToVertexLabels(const GRAPH& graph, const ELA& edge_labels, VLA& vertex_labels)
{
    struct mask
    {
        mask(const ELA& edge_labels) : edge_labels_(edge_labels)
        {}
        bool vertex(size_t i) const
            { return true; }
        bool edge(size_t i) const
            { return !edge_labels_[i]; }

        const ELA& edge_labels_;
    };

    andres::graph::DepthFirstSearchData<> dfs_data(graph.numberOfVertices());

    for (size_t i = 0, label = 0; i < graph.numberOfVertices(); ++i)
        if (!dfs_data.visited(i))
        {
            depthFirstSearch(
                graph,
                mask(edge_labels),
                i,
                [&](size_t v, bool& proceed, bool& add_neighbors)
                {
                    vertex_labels[v] = label;
                    proceed = true;
                    add_neighbors = true;
                },
                dfs_data);

            ++label;
        }
}

template<typename GRAPH>
void enforce_partial_solution(andres::View<size_t> const& partial_solution, nl_lmp::Problem<GRAPH>& problem, nl_lmp::Solution& solution)
{
    size_t max_cluster_index = 0;
    for (size_t i = 0; i < partial_solution.shape(0); ++i)
        if (partial_solution(i, 0) < problem.numberOfClasses())
        {
            for (size_t j = 0; j < problem.numberOfClasses(); ++j)
                problem.setUnaryCost(i, j, 1000.0);

            problem.setUnaryCost(i, partial_solution(i, 0), -1000.0);

            solution[i].classIndex = partial_solution(i, 0);

            solution[i].clusterIndex = partial_solution(i, 1);
            max_cluster_index = max(max_cluster_index, partial_solution(i, 1));
        }

    for (size_t i = 0; i < partial_solution.shape(0); ++i)
        if (partial_solution(i, 1) > max_cluster_index)
            solution[i].clusterIndex = max_cluster_index + 1;

    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
        for (size_t w = v + 1; w < problem.numberOfVertices(); ++w)
        {
            auto const v_class = partial_solution(v, 0);
            auto const w_class = partial_solution(w, 0);

            if (v_class < problem.numberOfClasses() && w_class < problem.numberOfClasses())
            {
                if (partial_solution(v, 1) != partial_solution(w, 1))
                    for (size_t k = 0; k < problem.numberOfClasses(); ++k)
                        for (size_t l = 0; l < problem.numberOfClasses(); ++l)
                        {
                            problem.setPairwiseCutCost(v, w, k, l, -1000.0);
                            problem.setPairwiseJoinCost(v, w, k, l, 1000.0);
                        }
                else
                    for (size_t k = 0; k < problem.numberOfClasses(); ++k)
                        for (size_t l = 0; l < problem.numberOfClasses(); ++l)
                        {
                            problem.setPairwiseCutCost(v, w, k, l, 1000.0);
                            problem.setPairwiseJoinCost(v, w, k, l, -1000.0);
                        }
            }
        }
}

template<typename GRAPH>
nl_lmp::Problem<GRAPH> make_problem(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, bool do_logit_transform, bool with_suppression = false)
{
    auto const numberOfVertices = unaries.shape(0);
    auto const numberOfClasses = unaries.shape(1);

    nl_lmp::Problem<GRAPH> problem(numberOfVertices, numberOfClasses + (with_suppression ? 1 : 0));
    
    andres::NegativeLogProbabilityRatio<double, double> logit(1e-7);

    for (size_t v = 0; v < numberOfVertices; ++v)
    {
        for (size_t c = 0; c < numberOfClasses; ++c)
        {
            auto p = unaries(v, c);

            if (!do_logit_transform)
                problem.setUnaryCost(v, c, p);
            else
            {
                if (p < 1e-7)
                    problem.setUnaryCost(v, c, 1000.0);
                else
                    problem.setUnaryCost(v, c, logit(p));
            }
        }
    }

    for (size_t e = 0; e < edges.shape(0); ++e)
        if (edges(e, 0) != edges(e, 1))
        {
            auto p = probabilities(e);

            problem.setPairwiseJoinCost(
                edges(e, 0),
                edges(e, 1),
                numberOfClasses == 1 ? 0 : edges(e, 2),
                numberOfClasses == 1 ? 0 : edges(e, 3),
                do_logit_transform ? logit(p) : p
            );
        }

    if (with_suppression)
        for (size_t e = 0; e < problem.liftedGraph().numberOfEdges(); ++e)
        {
            auto const v0 = problem.liftedGraph().vertexOfEdge(e, 0);
            auto const v1 = problem.liftedGraph().vertexOfEdge(e, 1);

            for (size_t k = 0; k < numberOfClasses; ++k)
            {
                problem.setPairwiseJoinCost(v0, v1, k, numberOfClasses, 1000.0, e);
                problem.setPairwiseJoinCost(v0, v1, numberOfClasses, k, 1000.0, e);
            }
        }

    cout << numberOfVertices << " vertices and " << numberOfClasses << " classes and " << problem.originalGraph().numberOfEdges() << " edges" << endl;

    return problem;
}

template<typename GRAPH>
andres::Marray<size_t> make_and_solve(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, andres::View<size_t> const& partial_solution, SolvingMethod method, bool do_logit_transform, bool with_suppression)
{
    auto problem = make_problem<GRAPH>(unaries, edges, probabilities, do_logit_transform, with_suppression);

    nl_lmp::Solution solution(problem.numberOfVertices());

    // initialize vertex labels with the most probable class
    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
    {
        double best_gain = problem.getUnaryCost(v, 0);
        size_t best_color = 0;

        for (size_t k = 1; k < problem.numberOfClasses() - (with_suppression ? 1 : 0); ++k)
            if (problem.getUnaryCost(v, k) < best_gain)
            {
                best_gain = problem.getUnaryCost(v, k);
                best_color = k;
            }

        solution[v].classIndex = best_color;
    }

    enforce_partial_solution(partial_solution, problem, solution);

    solution = nl_lmp::greedyAdditiveEdgeContraction(problem, solution);

    if (method == SolvingMethod::Simple)
        solution = nl_lmp::solve_alternating(problem, solution);
    else if (method == SolvingMethod::Joint)
        solution = nl_lmp::update_labels_and_multicut(problem, solution);
#ifdef WITH_GUROBI
    else if (method == SolvingMethod::ILP)
    {
        if (problem.numberOfClasses() > 1 || with_suppression)
            throw runtime_error("ILP solver is available only for the pure multicut problem");

        vector<double> edge_weights(problem.originalGraph().numberOfEdges());
        for (size_t e = 0; e < problem.originalGraph().numberOfEdges(); ++e)
        {
            auto const v0 = problem.originalGraph().vertexOfEdge(e, 0);
            auto const v1 = problem.originalGraph().vertexOfEdge(e, 1);

            edge_weights[e] = problem.getPairwiseCutCost(v0, v1, solution[v0].classIndex, solution[v1].classIndex, e) - problem.getPairwiseJoinCost(v0, v1, solution[v0].classIndex, solution[v1].classIndex, e);
        }

        vector<char> edge_labels(problem.originalGraph().numberOfEdges());

        andres::graph::multicut::ilp_callback<andres::ilp::GurobiCallback>(problem.originalGraph(), edge_weights, edge_labels, edge_labels);

        vector<size_t> vertex_labels(problem.numberOfVertices());

        edgeToVertexLabels(problem.originalGraph(), edge_labels, vertex_labels);

        for (size_t i = 0; i < vertex_labels.size(); ++i)
            solution[i].clusterIndex = vertex_labels[i];
    }
#endif

    double energy = .0;

    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
        energy += problem.getUnaryCost(v, solution[v].classIndex);

    for (size_t e = 0; e < problem.liftedGraph().numberOfEdges(); ++e)
    {
        auto const v0 = problem.liftedGraph().vertexOfEdge(e, 0);
        auto const v1 = problem.liftedGraph().vertexOfEdge(e, 1);

        if (solution[v0].clusterIndex == solution[v1].clusterIndex)
            energy += problem.getPairwiseJoinCost(v0, v1, solution[v0].classIndex, solution[v1].classIndex);
        else
            energy += problem.getPairwiseCutCost(v0, v1, solution[v0].classIndex, solution[v1].classIndex);
    }

    cout << "Objective: " << energy << endl;

    if (with_suppression)
        for (size_t v = 0; v < problem.numberOfVertices(); ++v)
            if (solution[v].classIndex == problem.numberOfClasses() - 1)
                solution[v].classIndex = numeric_limits<size_t>::max();

    andres::Marray<size_t> m({ solution.size(), 2 });
    
    for(size_t v = 0; v < solution.size(); ++v)
    {
        m(v, 0) = solution[v].classIndex;
        m(v, 1) = solution[v].clusterIndex;
    }

    return m;
}

andres::Marray<size_t> solve_nl_lmp_complete_graph(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, SolvingMethod method, bool do_logit_transform, bool with_suppression)
{
    andres::Marray<size_t> partial_solution( { unaries.shape(0), 2}, numeric_limits<size_t>::max() );

    return make_and_solve<andres::graph::CompleteGraph<>>(unaries, edges, probabilities, partial_solution, method, do_logit_transform, with_suppression);
}

andres::Marray<size_t> solve_nl_lmp_sparse_graph(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, SolvingMethod method, bool do_logit_transform, bool with_suppression)
{
    andres::Marray<size_t> partial_solution( { unaries.shape(0), 2}, numeric_limits<size_t>::max() );

    return make_and_solve<andres::graph::Graph<>>(unaries, edges, probabilities, partial_solution, method, do_logit_transform, with_suppression);
}

andres::Marray<size_t> solve_nl_lmp_complete_graph(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, andres::View<size_t> const& partial_solution, SolvingMethod method, bool do_logit_transform, bool with_suppression)
{
    return make_and_solve<andres::graph::CompleteGraph<>>(unaries, edges, probabilities, partial_solution, method, do_logit_transform, with_suppression);
}

andres::Marray<size_t> solve_nl_lmp_sparse_graph(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, andres::View<size_t> const& partial_solution, SolvingMethod method, bool do_logit_transform, bool with_suppression)
{
    return make_and_solve<andres::graph::Graph<>>(unaries, edges, probabilities, partial_solution, method, do_logit_transform, with_suppression);
}

double compute_objective(andres::View<double> const& unaries, andres::View<uint16_t> const& edges, andres::View<double> const& probabilities, andres::View<size_t> const& solution, bool do_logit_transform)
{
    auto const numberOfClasses = unaries.shape(1);

    auto problem = make_problem<andres::graph::Graph<>>(unaries, edges, probabilities, do_logit_transform);

    double energy = .0;

    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
        if (solution(v, 0) < numberOfClasses)
            energy += problem.getUnaryCost(v, solution(v, 0));

    for (size_t e = 0; e < problem.liftedGraph().numberOfEdges(); ++e)
    {
        auto const v0 = problem.liftedGraph().vertexOfEdge(e, 0);
        auto const v1 = problem.liftedGraph().vertexOfEdge(e, 1);

        if (solution(v0, 0) >= numberOfClasses || solution(v1, 0) >= numberOfClasses)
            continue;

        if (solution(v0, 1) == solution(v1, 1))
            energy += problem.getPairwiseJoinCost(v0, v1, solution(v0, 0), solution(v1, 0));
        else
            energy += problem.getPairwiseCutCost(v0, v1, solution(v0, 0), solution(v1, 0));
    }

    return energy;
}
