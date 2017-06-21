#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_KERNIGHAN_LIN_HXX
#define ANDRES_GRAPH_MULTICUT_KERNIGHAN_LIN_HXX

#include <iomanip>
#include <stack>
#include <stdexcept>
#include <set>
#include <vector>

#include "../complete-graph.hxx"



namespace andres {
namespace graph {
namespace multicut {

struct KernighanLinSettings
{
    size_t numberOfInnerIterations { std::numeric_limits<size_t>::max() };
    size_t numberOfOuterIterations { 100 };
    double epsilon { 1e-6 };
    bool verbose { false };
};

template<typename GRAPH, typename ECA, typename VERTEXLABELS>
inline
auto kernighanLin(const GRAPH& graph, const ECA& edge_costs, const VERTEXLABELS& input_vertex_labels, const KernighanLinSettings settings = KernighanLinSettings()) -> VERTEXLABELS
{
    struct Visitor
    {
        bool operator()(VERTEXLABELS const& vertex_labels) const
        {
            return true;
        }
    } visitor;

    return kernighanLin(graph, edge_costs, input_vertex_labels, visitor, settings);
}

template<typename GRAPH, typename ECA, typename VERTEXLABELS, typename VIS>
inline
auto kernighanLin(const GRAPH& graph, const ECA& edge_costs, const VERTEXLABELS& input_vertex_labels, VIS& visitor, const KernighanLinSettings settings = KernighanLinSettings()) -> VERTEXLABELS
{
    struct TwoCutBuffers
    {
        TwoCutBuffers(const GRAPH& graph) :
            differences(graph.numberOfVertices()),
            is_moved(graph.numberOfVertices()),
            referenced_by(graph.numberOfVertices()),
            vertex_labels(graph.numberOfVertices())
        {}

        std::vector<size_t> border;
        std::vector<double> differences;
        std::vector<char> is_moved;
        size_t max_not_used_label;
        std::vector<size_t> referenced_by;
        VERTEXLABELS vertex_labels;
    } buffer(graph);

    auto update_bipartition = [&](std::vector<size_t>& A, std::vector<size_t>& B)
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

                for (auto it = graph.adjacenciesFromVertexBegin(A[i]); it != graph.adjacenciesFromVertexEnd(A[i]); ++it)
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

        buffer.border.clear();
        
        for (auto a : A)
            if (buffer.referenced_by[a] > 0)
                buffer.border.push_back(a);

        for (auto b : B)
            if (buffer.referenced_by[b] > 0)
                buffer.border.push_back(b);

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
                size_t size = buffer.border.size();
                
                for (size_t i = 0; i < size; )
                    if (buffer.referenced_by[buffer.border[i]] == 0)
                        std::swap(buffer.border[i], buffer.border[--size]);
                    else
                    {
                        if (buffer.differences[buffer.border[i]] > m.difference)
                        {
                            m.v = buffer.border[i];
                            m.difference = buffer.differences[m.v];
                        }
                        
                        ++i;
                    }

                buffer.border.erase(buffer.border.begin() + size, buffer.border.end());
            }
            
            if (m.v == -1)
                break;

            const auto old_label = buffer.vertex_labels[m.v];

            if (old_label == label_A)
                m.new_label = label_B;
            else
                m.new_label = label_A;

            // update differences and references
            for (auto it = graph.adjacenciesFromVertexBegin(m.v); it != graph.adjacenciesFromVertexEnd(m.v); ++it)
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
                        buffer.border.push_back(it->vertex());
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

    double starting_energy = .0;

    // check if the input multicut labeling is valid
    for(size_t edge = 0; edge < graph.numberOfEdges(); ++edge)
    {
        auto v0 = graph.vertexOfEdge(edge, 0);
        auto v1 = graph.vertexOfEdge(edge, 1);

        if (input_vertex_labels[v0] != input_vertex_labels[v1])
            starting_energy += edge_costs[edge];
    }

    auto numberOfComponents = *std::max_element(input_vertex_labels.begin(), input_vertex_labels.end()) + 1;

    // build partitions
    std::vector<std::vector<size_t>> partitions(numberOfComponents);

    for (size_t i = 0; i < input_vertex_labels.size(); ++i)
    {
        partitions[input_vertex_labels[i]].push_back(i);
        buffer.vertex_labels[i] = input_vertex_labels[i];
    }
    
    buffer.max_not_used_label = partitions.size();

    if (settings.verbose)
    {
        std::cout << "Starting energy: " << starting_energy << std::endl;
        std::cout << std::setw(4) << "Iter" << std::setw(16) << "Total decrease" << std::setw(15) << "Pair updates" << std::setw(15) << "New sets" << std::setw(15) << "Num. of sets\n";
    }

    auto last_good_vertex_labels = buffer.vertex_labels;

    // auxillary array for BFS/DFS
    std::vector<char> visited(graph.numberOfVertices());

    // 1 if i-th partitioned changed since last iteration, 0 otherwise
    std::vector<char> changed(numberOfComponents, 1);

    // interatively update bipartition in order to minimize the total cost of the multicut
    for (size_t k = 0; k < settings.numberOfOuterIterations; ++k)
    {
        auto energy_decrease = .0;

        std::vector<std::set<size_t>> edges(numberOfComponents);
        for (size_t e = 0; e < graph.numberOfEdges(); ++e)
        {
            auto const v0 = buffer.vertex_labels[graph.vertexOfEdge(e, 0)];
            auto const v1 = buffer.vertex_labels[graph.vertexOfEdge(e, 1)];

            if (v0 != v1)
                edges[std::min(v0, v1)].insert(std::max(v0, v1));
        }

        for (size_t i = 0; i < numberOfComponents; ++i)
            if (!partitions[i].empty())
                for (auto j : edges[i])
                    if (!partitions[j].empty() && (changed[j] || changed[i]))
                    {
                        auto ret = update_bipartition(partitions[i], partitions[j]);

                        if (ret > settings.epsilon)
                            changed[i] = changed[j] = 1;

                        energy_decrease += ret;

                        if (partitions[i].size() == 0)
                            break;
                    }
        
        auto ee = energy_decrease;

        // remove partitions that became empty after the previous step
        auto new_end = std::partition(partitions.begin(), partitions.end(), [](const std::vector<size_t>& s) { return !s.empty(); });
        partitions.resize(new_end - partitions.begin());

        // try to intoduce new partitions
        for (size_t i = 0, p_size = partitions.size(); i < p_size; ++i)
        {
            if (!changed[i])
                continue;

            while (1)
            {
                std::vector<size_t> new_set;
                energy_decrease += update_bipartition(partitions[i], new_set);

                if (new_set.empty())
                    break;

                partitions.emplace_back(std::move(new_set));
            }
        }

        if (!visitor(buffer.vertex_labels))
            break;

        if (energy_decrease == .0)
            break;

        std::stack<size_t> S;
        
        std::fill(visited.begin(), visited.end(), 0);

        partitions.clear();
        numberOfComponents = 0;

        // do connected component labeling on the original graph and form new pфrtitions
        for (size_t i = 0; i < graph.numberOfVertices(); ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto label = buffer.vertex_labels[i];

                buffer.referenced_by[i] = numberOfComponents;

                partitions.emplace_back(std::vector<size_t>());
                partitions.back().push_back(i);

                while (!S.empty())
                {
                    auto v = S.top();
                    S.pop();

                    for (auto it = graph.adjacenciesFromVertexBegin(v); it != graph.adjacenciesFromVertexEnd(v); ++it)
                        if (buffer.vertex_labels[it->vertex()] == label && !visited[it->vertex()])
                        {
                            S.push(it->vertex());
                            visited[it->vertex()] = 1;
                            buffer.referenced_by[it->vertex()] = numberOfComponents;
                            partitions.back().push_back(it->vertex());
                        }
                }

                ++numberOfComponents;
            }

        buffer.vertex_labels = buffer.referenced_by;

        buffer.max_not_used_label = numberOfComponents;

        bool didnt_change = true;
        for (size_t i = 0; i < graph.numberOfEdges(); ++i)
        {
            auto const v0 = graph.vertexOfEdge(i, 0);
            auto const v1 = graph.vertexOfEdge(i, 1);

            auto edge_label = buffer.vertex_labels[v0] == buffer.vertex_labels[v1] ? 0 : 1;

            if (static_cast<bool>(edge_label) != (last_good_vertex_labels[v0] != last_good_vertex_labels[v1]))
                didnt_change = false;
        }

        if (didnt_change)
            break;

        // check if the shape of some partitions didn't change
        changed.resize(numberOfComponents);
        std::fill(changed.begin(), changed.end(), 0);

        std::fill(visited.begin(), visited.end(), 0);

        for (size_t i = 0; i < graph.numberOfVertices(); ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto label_new = buffer.vertex_labels[i];
                auto label_old = last_good_vertex_labels[i];

                while (!S.empty())
                {
                    auto v = S.top();
                    S.pop();

                    for (auto w = graph.verticesFromVertexBegin(v); w != graph.verticesFromVertexEnd(v); ++w)
                    {
                        if (last_good_vertex_labels[*w] == label_old && buffer.vertex_labels[*w] != label_new)
                            changed[label_new] = 1;

                        if (visited[*w])
                            continue;

                        if (buffer.vertex_labels[*w] == label_new)
                        {
                            S.push(*w);
                            visited[*w] = 1;

                            if (last_good_vertex_labels[*w] != label_old)
                                changed[label_new] = 1;
                        }
                    }
                }
            }

        last_good_vertex_labels = buffer.vertex_labels;

        if (settings.verbose)
            std::cout << std::setw(4) << k+1 << std::setw(16) << energy_decrease << std::setw(15) << ee << std::setw(15) << (energy_decrease - ee) << std::setw(14) << partitions.size() << std::endl;
    }

    return last_good_vertex_labels;
}

template<typename GraphVisitor, typename ECA, typename VERTEXLABELS>
inline
auto kernighanLin(const CompleteGraph<GraphVisitor>& graph, const ECA& edge_costs, const VERTEXLABELS& input_vertex_labels, const KernighanLinSettings settings = KernighanLinSettings()) -> VERTEXLABELS
{
    struct Visitor
    {
        bool operator()(VERTEXLABELS const& edge_labels) const
        {
            return true;
        }
    } visitor;

    return kernighanLin(graph, edge_costs, input_vertex_labels, visitor, settings);
}

template<typename GraphVisitor, typename ECA, typename VERTEXLABELS, typename VIS>
inline
auto kernighanLin(const CompleteGraph<GraphVisitor>& graph, const ECA& edge_costs, const VERTEXLABELS& input_vertex_labels, VIS& visitor, const KernighanLinSettings settings = KernighanLinSettings()) -> VERTEXLABELS
{
    struct Buffers
    {
        Buffers(const CompleteGraph<GraphVisitor>& graph) :
            differences(graph.numberOfVertices()),
            is_moved(graph.numberOfVertices())
        {}

        std::vector<double> differences;
        std::vector<char> is_moved;
    } buffer(graph);

    auto update_bipartition = [&](std::vector<size_t>& A, std::vector<size_t>& B)
    {
        if (A.empty())
            return .0;

        auto gain_from_merging = .0;

        // compute differences for set A
        for (size_t i = 0; i < A.size(); ++i)
        {
            double diff = .0;
            buffer.is_moved[A[i]] = 0;
            
            for (auto v : A)
                if (A[i] != v)
                    diff -= edge_costs[graph.findEdge(A[i], v).second];

            for (auto v : B)
            {
                diff += edge_costs[graph.findEdge(A[i], v).second];
                gain_from_merging += edge_costs[graph.findEdge(A[i], v).second];
            }

            buffer.differences[A[i]] = diff;
        }

        // compute differences for set B
        for (size_t i = 0; i < B.size(); ++i)
        {
            double diff = .0;
            buffer.is_moved[B[i]] = 0;

            for (auto v : B)
                if (B[i] != v)
                    diff -= edge_costs[graph.findEdge(B[i], v).second];

            for (auto v : A)
                diff += edge_costs[graph.findEdge(B[i], v).second];

            buffer.differences[B[i]] = diff;
        }

        struct Move
        {
            int v { -1 };
            double difference { std::numeric_limits<double>::lowest() };
            char new_label;
        };

        double cumulative_diff = .0;
        std::pair<double, size_t> max_move { std::numeric_limits<double>::lowest(), 0 };
        std::vector<Move> moves;
        
        for (size_t k = 0; k < settings.numberOfInnerIterations; ++k)
        {
            Move m;

            for (auto a : A)
                if (!buffer.is_moved[a] && buffer.differences[a] > m.difference)
                {
                    m.v = a;
                    m.difference = buffer.differences[a];
                    m.new_label = 'B';
                }

            for (auto b : B)
                if (!buffer.is_moved[b] && buffer.differences[b] > m.difference)
                {
                    m.v = b;
                    m.difference = buffer.differences[b];
                    m.new_label = 'A';
                }

            if (m.v == -1)
                break;

            // update differences
            if (m.new_label == 'B')
            {
                for (auto v : A)
                    if (v != m.v && !buffer.is_moved[v] && buffer.differences[v] > std::numeric_limits<double>::lowest())
                        buffer.differences[v] += 2.0*edge_costs[graph.findEdge(v, m.v).second];

                for (auto v : B)
                    if (!buffer.is_moved[v] && buffer.differences[v] > std::numeric_limits<double>::lowest())
                        buffer.differences[v] -= 2.0*edge_costs[graph.findEdge(v, m.v).second];
                
            }
            else
            {
                for (auto v : A)
                    if (!buffer.is_moved[v] && buffer.differences[v] > std::numeric_limits<double>::lowest())
                        buffer.differences[v] -= 2.0*edge_costs[graph.findEdge(v, m.v).second];

                for (auto v : B)
                    if (v != m.v && !buffer.is_moved[v] && buffer.differences[v] > std::numeric_limits<double>::lowest())
                        buffer.differences[v] += 2.0*edge_costs[graph.findEdge(v, m.v).second];
            }

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

            B.clear();

            return gain_from_merging;
        }
        if (max_move.first > settings.epsilon)
        {
            for (size_t i = max_move.second; i < moves.size(); ++i)
                buffer.is_moved[moves[i].v] = 0;

            A.erase(std::partition(A.begin(), A.end(), [&](size_t a) { return !buffer.is_moved[a]; }), A.end());
            B.erase(std::partition(B.begin(), B.end(), [&](size_t b) { return !buffer.is_moved[b]; }), B.end());

            for (size_t i = 0; i < max_move.second; ++i)
                // move vertex to the other set
                if (moves[i].new_label == 'B')
                    B.push_back(moves[i].v);
                else
                    A.push_back(moves[i].v);

            return max_move.first;
        }

        return .0;
    };

    double starting_energy = .0;

    // check if the input multicut labeling is valid
    for(size_t edge = 0; edge < graph.numberOfEdges(); ++edge)
    {
        auto v0 = graph.vertexOfEdge(edge, 0);
        auto v1 = graph.vertexOfEdge(edge, 1);

        if (input_vertex_labels[v0] != input_vertex_labels[v1])
            starting_energy += edge_costs[edge];
    }

    auto numberOfComponents = *std::max_element(input_vertex_labels.begin(), input_vertex_labels.end()) + 1;

    // build partitions
    std::vector<std::vector<size_t>> partitions(numberOfComponents);
    for (size_t i = 0; i < input_vertex_labels.size(); ++i)
        partitions[input_vertex_labels[i]].push_back(i);

    if (settings.verbose)
    {
        std::cout << "Starting energy: " << starting_energy << std::endl;
        std::cout << std::setw(4) << "Iter" << std::setw(16) << "Total decrease" << std::setw(15) << "Pair updates" << std::setw(15) << "New sets" << std::setw(15) << "Num. of sets\n";
    }

    // interatively update bipartition in order to minimize the total cost of the multicut
    for (size_t k = 0; k < settings.numberOfOuterIterations; ++k)
    {
        auto energy_decrease = .0;

        // update pairs of partitions
        for (size_t i = 0; i < partitions.size() - 1; ++i)
            for (auto j = i + 1; j < partitions.size(); ++j)
                if (!partitions[j].empty())
                    energy_decrease += update_bipartition(partitions[i], partitions[j]);

        // remove partitions that became empty after the previous step
        auto new_end = std::partition(partitions.begin(), partitions.end(), [](const std::vector<size_t>& s) { return !s.empty(); });
        partitions.resize(new_end - partitions.begin());

        auto ee = energy_decrease;

        // try to intoduce new partitions
        for (size_t i = 0, p_size = partitions.size(); i < p_size; ++i)       
            while (1)
            {
                std::vector<size_t> new_set;
                energy_decrease += update_bipartition(partitions[i], new_set);

                if (!new_set.empty())
                    partitions.emplace_back(std::move(new_set));
                else
                    break;
            }

        if (energy_decrease == .0)
            break;

        if (settings.verbose)
            std::cout << std::setw(4) << k+1 << std::setw(16) << energy_decrease << std::setw(15) << ee << std::setw(15) << (energy_decrease - ee) << std::setw(14) << partitions.size() << std::endl;
    }

    VERTEXLABELS vertex_labels(graph.numberOfVertices());
    for (size_t i = 0; i < partitions.size(); ++i)
        for (size_t j = 0; j < partitions[i].size(); ++j)
            vertex_labels[partitions[i][j]] = i;

    return vertex_labels;
}


// functions for back compatability with old interface that works with edge labels
template<typename GRAPH, typename ECA, typename ELA>
inline
void kernighanLin(const GRAPH& graph, const ECA& edge_costs, const ELA& input_edge_labels, ELA& output_edge_labels, const KernighanLinSettings settings = KernighanLinSettings())
{
    struct Visitor
    {
        bool operator()(std::vector<size_t> const& vertex_labels) const
        {
            return true;
        }
    } visitor;

    return kernighanLin(graph, edge_costs, input_edge_labels, output_edge_labels, visitor, settings);
}

template<typename GRAPH, typename ECA, typename ELA, typename VIS>
inline
void kernighanLin(const GRAPH& graph, const ECA& edge_costs, const ELA& input_edge_labels, ELA& output_edge_labels, VIS& visitor, const KernighanLinSettings settings = KernighanLinSettings())
{
    std::vector<size_t> vertex_labels(graph.numberOfVertices());
    std::vector<char> visited(graph.numberOfVertices());

    std::stack<size_t> S;
    for (size_t i = 0, label = 0; i < graph.numberOfVertices(); ++i)
        if (!visited[i])
        {
            S.push(i);
            visited[i] = 1;

            while (!S.empty())
            {
                auto v = S.top();
                S.pop();

                vertex_labels[v] = label;

                for (auto it = graph.adjacenciesFromVertexBegin(v); it != graph.adjacenciesFromVertexEnd(v); ++it)
                    if (!input_edge_labels[it->edge()] && !visited[it->vertex()])
                    {
                        S.push(it->vertex());
                        visited[it->vertex()] = 1;
                    }
            }

            ++label;
        }

    vertex_labels = kernighanLin(graph, edge_costs, vertex_labels, visitor, settings);

    for (size_t e = 0; e < graph.numberOfEdges(); ++e)
    {
        auto const v0 = graph.vertexOfEdge(e, 0);
        auto const v1 = graph.vertexOfEdge(e, 1);

        output_edge_labels[e] = (vertex_labels[v0] != vertex_labels[v1]) ? 1 : 0;
    }
}

template<typename GraphVisitor, typename ECA, typename ELA>
inline
void kernighanLin(const CompleteGraph<GraphVisitor>& graph, const ECA& edge_costs, const ELA& input_edge_labels, ELA& output_edge_labels, const KernighanLinSettings settings = KernighanLinSettings())
{
    struct Visitor
    {
        bool operator()(std::vector<size_t> const& vertex_labels) const
        {
            return true;
        }
    } visitor;

    return kernighanLin(graph, edge_costs, input_edge_labels, output_edge_labels, visitor, settings);
}

template<typename GraphVisitor, typename ECA, typename ELA, typename VIS>
inline
void kernighanLin(const CompleteGraph<GraphVisitor>& graph, const ECA& edge_costs, const ELA& input_edge_labels, ELA& output_edge_labels, VIS& visitor, const KernighanLinSettings settings = KernighanLinSettings())
{
    std::vector<size_t> vertex_labels(graph.numberOfVertices());
    std::vector<char> visited(graph.numberOfVertices());

    std::stack<size_t> S;
    for (size_t i = 0, label = 0; i < graph.numberOfVertices(); ++i)
        if (!visited[i])
        {
            S.push(i);
            visited[i] = 1;

            while (!S.empty())
            {
                auto v = S.top();
                S.pop();

                vertex_labels[v] = label;

                for (auto it = graph.adjacenciesFromVertexBegin(v); it != graph.adjacenciesFromVertexEnd(v); ++it)
                    if (!input_edge_labels[it->edge()] && !visited[it->vertex()])
                    {
                        S.push(it->vertex());
                        visited[it->vertex()] = 1;
                    }
            }

            ++label;
        }

    vertex_labels = kernighanLin(graph, edge_costs, vertex_labels, visitor, settings);

    for (size_t e = 0; e < graph.numberOfEdges(); ++e)
    {
        auto const v0 = graph.vertexOfEdge(e, 0);
        auto const v1 = graph.vertexOfEdge(e, 1);

        output_edge_labels[e] = (vertex_labels[v0] != vertex_labels[v1]) ? 1 : 0;
    }
}

} // of multicut
} // of graph
} // of andres

#endif
