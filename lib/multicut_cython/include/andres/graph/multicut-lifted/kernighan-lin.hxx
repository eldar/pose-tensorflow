#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_KERNIGHAN_LIN_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_KERNIGHAN_LIN_HXX

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <set>
#include <vector>
#include <stack>



namespace andres {
namespace graph {
namespace multicut_lifted {

struct KernighanLinSettings {
    size_t numberOfInnerIterations { std::numeric_limits<size_t>::max() };
    size_t numberOfOuterIterations { 100 };
    double epsilon { 1e-6 };
    bool verbose { false };
    bool introduce_new_sets { true };
};


template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename VERTEXLABELS>
auto kernighanLin(
    ORIGINAL_GRAPH const& original_graph,
    LIFTED_GRAPH const& lifted_graph,
    ECA const& edge_costs,
    VERTEXLABELS const& input_vertex_labels,
    KernighanLinSettings settings = KernighanLinSettings()) -> VERTEXLABELS
{
    struct Visitor {
        constexpr bool operator()(VERTEXLABELS const& vertex_labels)
            { return true; }
    } visitor;

    return kernighanLin(original_graph, lifted_graph, edge_costs, input_vertex_labels, visitor, settings);
}

template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename VERTEXLABELS, typename VIS>
inline
auto kernighanLin(
    ORIGINAL_GRAPH const& original_graph,
    LIFTED_GRAPH const& lifted_graph,
    ECA const& edge_costs,
    VERTEXLABELS const& input_vertex_labels,
    VIS& visitor,
    const KernighanLinSettings settings = KernighanLinSettings()) -> VERTEXLABELS
{
    struct Buffers {
        Buffers(const ORIGINAL_GRAPH& graph)
            :   differences(graph.numberOfVertices()),
                is_moved(graph.numberOfVertices()),
                referenced_by(graph.numberOfVertices()),
                vertex_labels(graph.numberOfVertices())
            {}
        std::vector<double> differences;
        std::vector<char> is_moved;
        size_t max_not_used_label;
        std::vector<size_t> referenced_by;
        VERTEXLABELS vertex_labels;
    } buffer(original_graph);

    auto update_bipartition = [&buffer, &original_graph, &lifted_graph, &edge_costs, &settings] (std::vector<size_t>& A, std::vector<size_t>& B)
    {
        struct Move
        {
            int v { -1 };
            double difference { std::numeric_limits<double>::lowest() };
            size_t new_label;
        };

        auto gain_from_merging = .0;

        auto compute_differences = [&](const std::vector<size_t>& A, size_t label_A, size_t label_B)
        {
            for (long int i = 0; i < A.size(); ++i)
            {
                double diffExt = .0;
                double diffInt = .0;
                size_t ref_cnt = 0;

                if (lifted_graph.numberOfEdges() != original_graph.numberOfEdges())
                {
                    for (auto it = lifted_graph.adjacenciesFromVertexBegin(A[i]); it != lifted_graph.adjacenciesFromVertexEnd(A[i]); ++it)
                    {
                        const auto lbl = buffer.vertex_labels[it->vertex()];

                        if (lbl == label_A)
                            diffInt += edge_costs[it->edge()];
                        else if (lbl == label_B)
                            diffExt += edge_costs[it->edge()];
                    }

                    for (auto it = original_graph.adjacenciesFromVertexBegin(A[i]); it != original_graph.adjacenciesFromVertexEnd(A[i]); ++it)
                        if (buffer.vertex_labels[it->vertex()] == label_B)
                            ++ref_cnt;
                }
                else
                    for (auto it = original_graph.adjacenciesFromVertexBegin(A[i]); it != original_graph.adjacenciesFromVertexEnd(A[i]); ++it)
                    {
                        const auto lbl = buffer.vertex_labels[it->vertex()];

                        if (lbl == label_A)
                            diffInt += edge_costs[it->edge()];
                        else if (lbl == label_B)
                        {
                            diffExt += edge_costs[it->edge()];

                            ++ref_cnt;
                        }
                    }

                buffer.differences[A[i]] = diffExt - diffInt;
                buffer.referenced_by[A[i]] = ref_cnt;
                buffer.is_moved[A[i]] = 0;

                gain_from_merging += diffExt;
            }
        };


        if (A.empty())
            return .0;
        
        auto label_A = buffer.vertex_labels[A[0]];
        auto label_B = (!B.empty()) ? buffer.vertex_labels[B[0]] : buffer.max_not_used_label;

        compute_differences(A, label_A, label_B);
        compute_differences(B, label_B, label_A);

        gain_from_merging /= 2.0;

        std::vector<size_t> border;

        for (auto a : A)
            if (buffer.referenced_by[a] > 0)
                border.push_back(a);

        for (auto b : B)
            if (buffer.referenced_by[b] > 0)
                border.push_back(b);

        std::vector<Move> moves;
        double cumulative_diff = .0;
        std::pair<double, size_t> max_move { std::numeric_limits<double>::lowest(), 0 };

        for (size_t k = 0; k < settings.numberOfInnerIterations; ++k)
        {
            Move m;

            if (B.empty() && k == 0)
            {
                for (auto a : A)
                    if (buffer.differences[a] > m.difference)
                    {
                        m.v = a;
                        m.difference = buffer.differences[a];
                    }
            }
            else
            {
                size_t size = border.size();
                
                for (size_t i = 0; i < size; )
                    if (buffer.referenced_by[border[i]] == 0)
                        std::swap(border[i], border[--size]);
                    else
                    {
                        if (buffer.differences[border[i]] > m.difference)
                        {
                            m.v = border[i];
                            m.difference = buffer.differences[m.v];
                        }
                        
                        ++i;
                    }

                border.erase(border.begin() + size, border.end());
            }

            if (m.v == -1)
                break;

            const auto old_label = buffer.vertex_labels[m.v];

            if (old_label == label_A)
                m.new_label = label_B;
            else
                m.new_label = label_A;

            // update differences and references
            if (lifted_graph.numberOfEdges() != original_graph.numberOfEdges())
            {
                for (auto it = lifted_graph.adjacenciesFromVertexBegin(m.v); it != lifted_graph.adjacenciesFromVertexEnd(m.v); ++it)
                {
                    if (buffer.is_moved[it->vertex()])
                        continue;

                    const auto lbl = buffer.vertex_labels[it->vertex()];
                    
                    // edge to an element of the new set
                    if (lbl == m.new_label)
                        buffer.differences[it->vertex()] -= 2.0*edge_costs[it->edge()];
                    // edge to an element of the old set
                    else if (lbl == old_label)
                        buffer.differences[it->vertex()] += 2.0*edge_costs[it->edge()];
                }

                for (auto it = original_graph.adjacenciesFromVertexBegin(m.v); it != original_graph.adjacenciesFromVertexEnd(m.v); ++it)
                {
                    if (buffer.is_moved[it->vertex()])
                        continue;

                    const auto lbl = buffer.vertex_labels[it->vertex()];

                    // edge to an element of the new set
                    if (lbl == m.new_label)
                        --buffer.referenced_by[it->vertex()];
                    // edge to an element of the old set
                    else if (lbl == old_label)
                    {
                        ++buffer.referenced_by[it->vertex()];

                        if (buffer.referenced_by[it->vertex()] == 1)
                            border.push_back(it->vertex());
                    }
                }
            }
            else
                for (auto it = original_graph.adjacenciesFromVertexBegin(m.v); it != original_graph.adjacenciesFromVertexEnd(m.v); ++it)
                {
                    if (buffer.is_moved[it->vertex()])
                        continue;

                    const auto lbl = buffer.vertex_labels[it->vertex()];
                    
                    // edge to an element of the new set
                    if (lbl == m.new_label)
                    {
                        buffer.differences[it->vertex()] -= 2.0*edge_costs[it->edge()];

                        --buffer.referenced_by[it->vertex()];
                    }
                    // edge to an element of the old set
                    else if (lbl == old_label)
                    {
                        buffer.differences[it->vertex()] += 2.0*edge_costs[it->edge()];

                        ++buffer.referenced_by[it->vertex()];

                        if (buffer.referenced_by[it->vertex()] == 1)
                            border.push_back(it->vertex());
                    }
                }

            buffer.vertex_labels[m.v] = m.new_label;
            buffer.referenced_by[m.v] = 0;
            buffer.differences[m.v] = std::numeric_limits<double>::lowest();
            buffer.is_moved[m.v] = 1;
            moves.push_back(m);

            cumulative_diff += m.difference;

            if (cumulative_diff > max_move.first)
                max_move = std::make_pair(cumulative_diff, moves.size());
        }

        
        if (gain_from_merging > max_move.first && gain_from_merging > settings.epsilon)
        {
            A.insert(A.end(), B.begin(), B.end());

            for (auto a : A)
                buffer.vertex_labels[a] = label_A;

            for (auto b : B)
                buffer.vertex_labels[b] = label_A;

            B.clear();

            return gain_from_merging;
        }
        else if (max_move.first > settings.epsilon)
        {
            // revert some changes
            for (size_t i = max_move.second; i < moves.size(); ++i)
            {
                buffer.is_moved[moves[i].v] = 0;

                if (moves[i].new_label == label_B)
                    buffer.vertex_labels[moves[i].v] = label_A;
                else
                    buffer.vertex_labels[moves[i].v] = label_B;
            }

            // make sure that this is unique label
            if (B.empty())
                ++buffer.max_not_used_label;

            A.erase(std::partition(A.begin(), A.end(), [&](size_t a) { return !buffer.is_moved[a]; }), A.end());
            B.erase(std::partition(B.begin(), B.end(), [&](size_t b) { return !buffer.is_moved[b]; }), B.end());

            for (size_t i = 0; i < max_move.second; ++i)
                // move vertex to the other set
                if (moves[i].new_label == label_B)
                    B.push_back(moves[i].v);
                else
                    A.push_back(moves[i].v);

            return max_move.first;
        }
        else
            for (size_t i = 0; i < moves.size(); ++i)
                if (moves[i].new_label == label_B)
                    buffer.vertex_labels[moves[i].v] = label_A;
                else
                    buffer.vertex_labels[moves[i].v] = label_B;

        return .0;
    };

    auto compute_obj_value = [&lifted_graph, &edge_costs] (std::vector<size_t> const& vertex_labels)
    {
        double obj_value = .0;

        for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
        {
            auto const v0 = lifted_graph.vertexOfEdge(i, 0);
            auto const v1 = lifted_graph.vertexOfEdge(i, 1);

            if (vertex_labels[v0] != vertex_labels[v1])
                obj_value += edge_costs[i];
        }

        return obj_value;
    };

    auto form_partitions = [] (std::vector<size_t> const& vertex_labels)
    {
        auto const number_of_components = *std::max_element(vertex_labels.begin(), vertex_labels.end()) + 1;

        std::vector<std::vector<size_t>> partitions(number_of_components);

        for (size_t i = 0; i < vertex_labels.size(); ++i)
            partitions[vertex_labels[i]].push_back(i);
        
        return partitions;
    };

    auto project_onto_feasable_set = [&original_graph] (std::vector<size_t> const& vertex_labels)
    {
        std::stack<size_t> S;
        std::vector<char> visited(original_graph.numberOfVertices());

        std::vector<size_t> feasable_vertex_labels(vertex_labels.size());

        // do connected component labeling on the original graph
        for (size_t i = 0, new_label = 0; i < original_graph.numberOfVertices(); ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto const label = vertex_labels[i];

                while (!S.empty())
                {
                    auto const v = S.top();
                    S.pop();

                    feasable_vertex_labels[v] = new_label;

                    for (auto it = original_graph.adjacenciesFromVertexBegin(v); it != original_graph.adjacenciesFromVertexEnd(v); ++it)
                        if (vertex_labels[it->vertex()] == label && !visited[it->vertex()])
                        {
                            S.push(it->vertex());
                            visited[it->vertex()] = 1;
                        }
                }

                ++new_label;
            }

        return feasable_vertex_labels;
    };

    auto mark_partitions_that_changed_shape = [&original_graph] (std::vector<size_t> const& previous_vertex_labels, std::vector<size_t> const& current_vertex_labels)
    {
        auto const number_of_components = *std::max_element(current_vertex_labels.begin(), current_vertex_labels.end()) + 1;

        std::vector<char> changed(number_of_components);

        std::stack<size_t> S;
        std::vector<char> visited(original_graph.numberOfVertices());
        
        for (size_t i = 0; i < original_graph.numberOfVertices(); ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto const label_new = current_vertex_labels[i];
                auto const label_old = previous_vertex_labels[i];

                while (!S.empty())
                {
                    auto const v = S.top();
                    S.pop();

                    for (auto w = original_graph.verticesFromVertexBegin(v); w != original_graph.verticesFromVertexEnd(v); ++w)
                    {
                        if (
                            previous_vertex_labels[*w] == label_old && current_vertex_labels[*w] != label_new
                            ||
                            previous_vertex_labels[*w] != label_old && current_vertex_labels[*w] == label_new
                            )
                            changed[label_new] = 1;

                        if (visited[*w])
                            continue;

                        if (current_vertex_labels[*w] == label_new)
                        {
                            S.push(*w);
                            visited[*w] = 1;

                            if (previous_vertex_labels[*w] != label_old)
                                changed[label_new] = 1;
                        }
                    }
                }
            }

        return changed;
    };



    buffer.vertex_labels = project_onto_feasable_set(input_vertex_labels);

    auto current_obj_value = compute_obj_value(buffer.vertex_labels);

    auto partitions = form_partitions(buffer.vertex_labels);
    buffer.max_not_used_label = partitions.size();

    if (settings.verbose)
    {
        std::cout << "Starting number of segments: " << partitions.size() << std::endl;
        std::cout << "Starting energy: " << std::fixed << std::setprecision(4) << current_obj_value << std::endl;
        std::cout << std::setw(4) << "Iter" << std::setw(16) << "Obj. value" << std::setw(15) << "Pair updates" << std::setw(15) << "New sets" << std::setw(10) << "# sets\n";
    }

    auto last_good_vertex_labels = buffer.vertex_labels;
    auto vertex_labels_for_changed_thing = buffer.vertex_labels;

    for (size_t k = 0; k < settings.numberOfOuterIterations; ++k)
    {
        // 1 if i-th partitioned changed since last iteration, 0 otherwise
        std::vector<char> changed(partitions.size(), 1);

        if (k > 0)
            changed = mark_partitions_that_changed_shape(vertex_labels_for_changed_thing, buffer.vertex_labels);

        std::vector<std::set<size_t>> edges(partitions.size());
        for (size_t e = 0; e < original_graph.numberOfEdges(); ++e)
        {
            auto const v0 = buffer.vertex_labels[original_graph.vertexOfEdge(e, 0)];
            auto const v1 = buffer.vertex_labels[original_graph.vertexOfEdge(e, 1)];

            if (v0 != v1)
                edges[std::min(v0, v1)].insert(std::max(v0, v1));
        }

        for (size_t i = 0; i < partitions.size(); ++i)
            if (!partitions[i].empty())
                for (auto j = edges[i].begin(); j != edges[i].end(); ++j)
                    if (!partitions[*j].empty() && (changed[*j] || changed[i]))
                    {
                        update_bipartition(partitions[i], partitions[*j]);

                        if (partitions[i].size() == 0)
                            break;
                    }


        buffer.vertex_labels = project_onto_feasable_set(buffer.vertex_labels);

        partitions = form_partitions(buffer.vertex_labels);
        buffer.max_not_used_label = partitions.size();

        auto pair_updates_decrease = current_obj_value - compute_obj_value(buffer.vertex_labels);
        current_obj_value -= pair_updates_decrease;



        if (k > 0)
            changed = mark_partitions_that_changed_shape(vertex_labels_for_changed_thing, buffer.vertex_labels);
        else
            changed.resize(partitions.size(), 1);

        // try to intoduce new partitions
        if (settings.introduce_new_sets)
            for (size_t i = 0; i < partitions.size(); ++i)
            {
                if (!changed[i])
                    continue;

                while (1)
                {
                    std::vector<size_t> new_set;
                    update_bipartition(partitions[i], new_set);
                    
                    if (new_set.empty())
                        break;
                }
            }


        buffer.vertex_labels = project_onto_feasable_set(buffer.vertex_labels);

        partitions = form_partitions(buffer.vertex_labels);
        buffer.max_not_used_label = partitions.size();

        auto new_sets_decrease = current_obj_value - compute_obj_value(buffer.vertex_labels);
        current_obj_value -= new_sets_decrease;


        if (new_sets_decrease + pair_updates_decrease < std::numeric_limits<double>::epsilon())
            break;

        bool didnt_change = true;
        for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
        {
            auto v0 = lifted_graph.vertexOfEdge(i, 0);
            auto v1 = lifted_graph.vertexOfEdge(i, 1);

            auto edge_label = buffer.vertex_labels[v0] == buffer.vertex_labels[v1] ? 0 : 1;

            if (static_cast<bool>(edge_label) != (last_good_vertex_labels[v0] != last_good_vertex_labels[v1]))
                didnt_change = false;
        }

        if (didnt_change)
            break;

        vertex_labels_for_changed_thing = last_good_vertex_labels;

        last_good_vertex_labels = buffer.vertex_labels;

        if (!visitor(last_good_vertex_labels))
            break;

        if (settings.verbose)
            std::cout << std::setw(4) << k+1 << std::setw(16) << current_obj_value << std::setw(15) << pair_updates_decrease << std::setw(15) << new_sets_decrease << std::setw(10) << partitions.size() << std::endl;
    }

    if (settings.verbose)
        std::cout << "Final objective: " << compute_obj_value(last_good_vertex_labels) << std::endl;

    return last_good_vertex_labels;
}



// functions for back compatability with old interface that works with edge labels
template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename ELA>
inline
void kernighanLin(const ORIGINAL_GRAPH& original_graph, const LIFTED_GRAPH& lifted_graph, const ECA& edge_costs, const ELA& input_edge_labels, ELA& output_edge_labels, const KernighanLinSettings settings = KernighanLinSettings())
{
    struct Visitor
    {
        bool operator()(std::vector<size_t> const& vertex_labels) const
        {
            return true;
        }
    } visitor;

    return kernighanLin(original_graph, lifted_graph, edge_costs, input_edge_labels, output_edge_labels, visitor, settings);
}

template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename ELA, typename VIS>
inline
void kernighanLin(const ORIGINAL_GRAPH& original_graph, const LIFTED_GRAPH& lifted_graph, const ECA& edge_costs, const ELA& input_edge_labels, ELA& output_edge_labels, VIS& visitor, const KernighanLinSettings settings = KernighanLinSettings())
{
    std::vector<size_t> vertex_labels(lifted_graph.numberOfVertices());
    std::vector<char> visited(lifted_graph.numberOfVertices());

    std::stack<size_t> S;
    for (size_t i = 0, label = 0; i < lifted_graph.numberOfVertices(); ++i)
        if (!visited[i])
        {
            S.push(i);
            visited[i] = 1;

            while (!S.empty())
            {
                auto v = S.top();
                S.pop();

                vertex_labels[v] = label;

                for (auto it = lifted_graph.adjacenciesFromVertexBegin(v); it != lifted_graph.adjacenciesFromVertexEnd(v); ++it)
                    if (!input_edge_labels[it->edge()] && !visited[it->vertex()])
                    {
                        S.push(it->vertex());
                        visited[it->vertex()] = 1;
                    }
            }

            ++label;
        }

    vertex_labels = kernighanLin(original_graph, lifted_graph, edge_costs, vertex_labels, visitor, settings);

    for (size_t e = 0; e < lifted_graph.numberOfEdges(); ++e)
    {
        auto const v0 = lifted_graph.vertexOfEdge(e, 0);
        auto const v1 = lifted_graph.vertexOfEdge(e, 1);

        output_edge_labels[e] = (vertex_labels[v0] != vertex_labels[v1]) ? 1 : 0;
    }
}

}
}
}

#endif
