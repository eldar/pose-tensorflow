#pragma once
#ifndef NL_LMP_SOLVE_JOINT_HXX
#define NL_LMP_SOLVE_JOINT_HXX

#include <algorithm>
#include <map>
#include <stack>
#include <unordered_set>
#include <vector>

#include "detail/compute-objective.hxx"
#include "detail/update-class-labels.hxx"



namespace nl_lmp
{

// gets a problem and some initial valid solution and solves the problem
// the assumption is that vertex cluster labels are distinct for distinct connected components
template<typename GRAPH>
inline
Solution update_labels_and_multicut(Problem<GRAPH> const& problem, Solution const& input)
{
    struct EmptyVisitor
    {
        void operator()(Solution const& solution) const
        {}

        void operator()(std::vector<size_t> const& cluster_labels, std::vector<size_t> const& class_labels) const
        {}
    } visitor;

    return update_labels_and_multicut(problem, input, visitor);
}

// gets a problem and some initial valid solution and solves the problem
// the assumption is that vertex cluster labels are distinct for distinct connected components
// this version also gets a visitor using which one can closely follow the optimization process
template<typename GRAPH, typename VISITOR>
inline
Solution update_labels_and_multicut(Problem<GRAPH> const& problem, Solution const& input, VISITOR& visitor)
{
    struct TwoCutBuffers
    {
        TwoCutBuffers(size_t n, size_t number_of_classes) :
            differences(n*number_of_classes),
            is_moved(n),
            referenced_by(n),
            cluster_labels(n),
            class_labels(n)
        {}

        // constains vertices of the boundary between 2 segments
        std::vector<size_t> border;
        // for each vertex contains number_of_classes possible gains
        std::vector<double> differences;
        std::vector<char> is_moved;
        std::vector<size_t> referenced_by;
        std::vector<size_t> cluster_labels;
        std::vector<size_t> class_labels;
        size_t max_not_used_label;
    } buffer(problem.numberOfVertices(), problem.numberOfClasses());

    // gets as input 2 arays containing vertices of 2 neighboring sets
    // updates the boundary between them and returns the expected gain
    auto update_bipartition = [&problem, &buffer] (std::vector<size_t>& A, std::vector<size_t>& B)
    {
        auto const number_of_classes = problem.numberOfClasses();

        // for each vertex in A computes the gain of moving (or not) this vertex into set B with, possibly, changing its label
        auto compute_differences = [&problem, &buffer, &number_of_classes](std::vector<size_t> const& A, size_t label_A, size_t label_B)
        {
            for (long int i = 0; i < A.size(); ++i)
            {
                auto const v = A[i];
                auto const v_label = buffer.class_labels[v];

                // for each 'v' there are 4 options:
                // 1) don't move, but change class label
                // 2) move, but do not change class label
                // 3) move and change class label
                // 4) do nothing
                
                // gains of changing class label
                // does not depend on the cut change
                for (size_t k = 0; k < number_of_classes; ++k)
                    buffer.differences[v*number_of_classes + k] = problem.getUnaryCost(v, v_label) - problem.getUnaryCost(v, k);
                
                // gains of changing cluster label, i.e. moving node from A to B
                // must be computed for all possible class labels of 'v', including no change
                for (auto it = problem.liftedGraph().adjacenciesFromVertexBegin(v); it != problem.liftedGraph().adjacenciesFromVertexEnd(v); ++it)
                {
                    auto const lbl = buffer.cluster_labels[it->vertex()];
                    auto const w = it->vertex();
                    auto const w_label = buffer.class_labels[w];

                    if (lbl == label_A)
                        // node in former set A
                        for (size_t k = 0; k < number_of_classes; ++k)
                            buffer.differences[v*number_of_classes + k] += problem.getPairwiseJoinCost(v, w, v_label, w_label, it->edge()) - problem.getPairwiseCutCost(v, w, k, w_label, it->edge());
                    else if (lbl == label_B)
                        // node in new set B
                        for (size_t k = 0; k < number_of_classes; ++k)
                            buffer.differences[v*number_of_classes + k] += problem.getPairwiseCutCost(v, w, v_label, w_label, it->edge()) - problem.getPairwiseJoinCost(v, w, k, w_label, it->edge());
                    else
                        // unlike in ordinary KL, due to possible changes of label of 'v', we must consider edges that go to neither A nor B
                        for (size_t k = 0; k < number_of_classes; ++k)
                            buffer.differences[v*number_of_classes + k] += problem.getPairwiseCutCost(v, w, v_label, w_label, it->edge()) - problem.getPairwiseCutCost(v, w, k, w_label, it->edge());
                }

                // this just finds out how many neighbors in B node 'v' has
                size_t ref_cnt = 0;
                for (auto it = problem.originalGraph().adjacenciesFromVertexBegin(v); it != problem.originalGraph().adjacenciesFromVertexEnd(v); ++it)
                    if (buffer.cluster_labels[it->vertex()] == label_B)
                        ++ref_cnt;

                buffer.referenced_by[v] = ref_cnt;
                buffer.is_moved[v] = 0;
            }
        };


        if (A.empty())
            return .0;

        auto label_A = buffer.cluster_labels[A[0]];
        // if B is empty, take the next free cluster label
        auto label_B = (!B.empty()) ? buffer.cluster_labels[B[0]] : buffer.max_not_used_label;
        
        compute_differences(A, label_A, label_B);
        compute_differences(B, label_B, label_A);

        double gain_from_merging = .0;
        for (auto a : A)
            if (buffer.referenced_by[a] > 0) 
            {
                auto const a_label = buffer.class_labels[a];

                for (auto it = problem.liftedGraph().adjacenciesFromVertexBegin(a); it != problem.liftedGraph().adjacenciesFromVertexEnd(a); ++it)
                {
                    auto const lbl = buffer.cluster_labels[it->vertex()];
                    auto const w = it->vertex();
                    auto const w_label = buffer.class_labels[w];

                    if (lbl == label_B)
                        gain_from_merging += problem.getPairwiseCutCost(a, w, a_label, w_label, it->edge()) - problem.getPairwiseJoinCost(a, w, a_label, w_label, it->edge());
                }
            }

        buffer.border.clear();
        
        // exactly maintain only the vertices in A and B that have neighbors in the other sets
        for (auto a : A)
            if (buffer.referenced_by[a] > 0)
                buffer.border.push_back(a);

        for (auto b : B)
            if (buffer.referenced_by[b] > 0)
                buffer.border.push_back(b);

        struct Move
        {
            int v { -1 };
            double difference { std::numeric_limits<double>::lowest() };
            size_t new_label;
            size_t new_class_index;
            size_t old_class_index;
        };

        std::vector<Move> moves;
        double cumulative_gain = .0;
        std::pair<double, size_t> max_move { std::numeric_limits<double>::lowest(), 0 };
        // assess moving all vertices of both sets
        // threshold on the best discovered cumulative gain
        for (size_t z = 0; ; ++z)
        {
            Move m;
            
            if (B.empty() && z == 0)
            {
                // if it's the very first iteration nad set B is empty, ther is no border yet, so examine all vertices in A
                for (auto a : A)
                    for (size_t k = 0; k < number_of_classes; ++k)
                        if (buffer.differences[a*number_of_classes + k] > m.difference)
                        {
                            m.v = a;
                            m.new_class_index = k;
                            m.difference = buffer.differences[a*number_of_classes + k];
                        }
            }
            else
            {
                size_t size = buffer.border.size();
                
                // this goes through all vertices that are on the boundary between A and B and chooses the next to move
                // at the same time it removes from 'border' vertices that stopped being in the boundary
                // this is very C++-specific, might not work for other languages
                for (size_t i = 0; i < size; )
                    if (buffer.referenced_by[buffer.border[i]] == 0)
                        // invariant assumption:
                        // if the vertex is not referenced by any other, then it's not on the boundary
                        // since boundary is not ordered we can just swap with the last element
                        std::swap(buffer.border[i], buffer.border[--size]);
                    else
                    {
                        // every vertex may also change it's class label, so examine all the class labels options
                        for (size_t k = 0; k < number_of_classes; ++k)
                            if (buffer.differences[buffer.border[i]*number_of_classes + k] > m.difference)
                            {
                                m.v = buffer.border[i];
                                m.new_class_index = k;
                                m.difference = buffer.differences[m.v*number_of_classes + k];
                            }
                        
                        ++i;
                    }
                
                // remove the vertices that stopped being on the boundary
                buffer.border.erase(buffer.border.begin() + size, buffer.border.end());
            }
            
            // nothing to move, break
            if (m.v == -1)
                break;
            
            // old cluster label
            auto const old_label = buffer.cluster_labels[m.v];
            
            // move the vertex to the other set based on it's current cluster label
            if (old_label == label_A)
                m.new_label = label_B;
            else
                m.new_label = label_A;

            m.old_class_index = buffer.class_labels[m.v];

            // update differences and references
            for (auto it = problem.liftedGraph().adjacenciesFromVertexBegin(m.v); it != problem.liftedGraph().adjacenciesFromVertexEnd(m.v); ++it)
            {
                // skip verticecs that have been moved
                if (buffer.is_moved[it->vertex()])
                    continue;

                auto const w = it->vertex();
                auto const w_label = buffer.class_labels[w];

                auto const lbl = buffer.cluster_labels[w];
                if (lbl == m.new_label)
                    // edge to an vertex of the new set
                    // update all its differences based on the new situation
                    for (size_t k = 0; k < number_of_classes; ++k)
                    {
                        // cancel old contribution
                        buffer.differences[w*number_of_classes + k] -= problem.getPairwiseCutCost(m.v, w, m.old_class_index, w_label, it->edge()) - problem.getPairwiseJoinCost(m.v, w, m.old_class_index, k, it->edge());

                        // add new contribution
                        buffer.differences[w*number_of_classes + k] += problem.getPairwiseJoinCost(m.v, w, m.new_class_index, w_label, it->edge()) - problem.getPairwiseCutCost(m.v, w, m.new_class_index, k, it->edge());
                    }
                else if (lbl == old_label)
                    // edge to an element of the old set
                    // update all its differences based on the new situation
                    for (size_t k = 0; k < number_of_classes; ++k)
                    {
                        // cancel old contribution
                        buffer.differences[w*number_of_classes + k] -= problem.getPairwiseJoinCost(m.v, w, m.old_class_index, w_label, it->edge()) - problem.getPairwiseCutCost(m.v, w, m.old_class_index, k, it->edge());

                        // add new contribution
                        buffer.differences[w*number_of_classes + k] += problem.getPairwiseCutCost(m.v, w, m.new_class_index, w_label, it->edge()) - problem.getPairwiseJoinCost(m.v, w, m.new_class_index, k, it->edge());
                    }
            }
            
            // this simply updates the boundary between two sets
            for (auto it = problem.originalGraph().adjacenciesFromVertexBegin(m.v); it != problem.originalGraph().adjacenciesFromVertexEnd(m.v); ++it)
            {
                if (buffer.is_moved[it->vertex()])
                    continue;

                const auto lbl = buffer.cluster_labels[it->vertex()];

                // edge to an element of the new set
                if (lbl == m.new_label)
                    // if the counter is 0, in the next iteration this vertex will be removed from the boundary
                    // the countaer is guaranteed to be >= 0
                    --buffer.referenced_by[it->vertex()];
                // edge to an element of the old set
                else if (lbl == old_label)
                {
                    ++buffer.referenced_by[it->vertex()];

                    if (buffer.referenced_by[it->vertex()] == 1)
                        // it it's the first time the vertex is referenced, add it to the boundary
                        buffer.border.push_back(it->vertex());
                }
            }

            buffer.cluster_labels[m.v] = m.new_label;
            buffer.referenced_by[m.v] = 0;
            buffer.is_moved[m.v] = 1;
            
            buffer.class_labels[m.v] = m.new_class_index;

            moves.push_back(m);

            cumulative_gain += m.difference;
            if (cumulative_gain > max_move.first)
                max_move = std::make_pair(cumulative_gain, moves.size());
        }

        if (gain_from_merging > max_move.first && gain_from_merging > 1e-6)
        {
            // revert all class label changes
            for (size_t i = 0; i < moves.size(); ++i)
                buffer.class_labels[moves[i].v] = moves[i].old_class_index;

            // if it's better to merge two partitions, just put everything into set A
            A.insert(A.end(), B.begin(), B.end());

            for (auto a : A)
                buffer.cluster_labels[a] = label_A;

            for (auto b : B)
                buffer.cluster_labels[b] = label_A;

            B.clear();

            return gain_from_merging;
        }
        else if (max_move.first > 1e-6)
        {
            // revert some changes
            for (size_t i = max_move.second; i < moves.size(); ++i)
            {
                buffer.is_moved[moves[i].v] = 0;
                
                // restore original cluster labels
                if (moves[i].new_label == label_B)
                    buffer.cluster_labels[moves[i].v] = label_A;
                else
                    buffer.cluster_labels[moves[i].v] = label_B;

                // restore original class labels
                buffer.class_labels[moves[i].v] = moves[i].old_class_index;
            }

            if (B.empty())
                ++buffer.max_not_used_label;
            
            // if 'v' is indeed moved, then buffer.is_moved[v] == 0, this is invariance
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
            // revert all changes
            for (size_t i = 0; i < moves.size(); ++i)
            {
                if (moves[i].new_label == label_B)
                    buffer.cluster_labels[moves[i].v] = label_A;
                else
                    buffer.cluster_labels[moves[i].v] = label_B;

                buffer.class_labels[moves[i].v] = moves[i].old_class_index;
            }

        return .0;
    };

    auto compute_obj_value = [&problem](std::vector<size_t> const& class_labels, std::vector<size_t> const& cluster_labels)
    {
        double obj_value = .0;

        for (size_t v = 0; v < problem.numberOfVertices(); ++v)
            obj_value += problem.getUnaryCost(v, class_labels[v]);

        for (size_t e = 0; e < problem.liftedGraph().numberOfEdges(); ++e)
        {
            auto const v0 = problem.liftedGraph().vertexOfEdge(e, 0);
            auto const v1 = problem.liftedGraph().vertexOfEdge(e, 1);

            if (cluster_labels[v0] != cluster_labels[v1])
                obj_value += problem.getPairwiseCutCost(v0, v1, class_labels[v0], class_labels[v1]);
            else
                obj_value += problem.getPairwiseJoinCost(v0, v1, class_labels[v0], class_labels[v1]);
        }

        return obj_value;
    };

    auto form_partitions = [] (std::vector<size_t> const& cluster_labels)
    {
        auto const number_of_components = *std::max_element(cluster_labels.begin(), cluster_labels.end()) + 1;

        std::vector<std::vector<size_t>> partitions(number_of_components);

        for (size_t i = 0; i < cluster_labels.size(); ++i)
            partitions[cluster_labels[i]].push_back(i);
        
        return partitions;
    };

    auto mark_partitions_that_changed = [&problem] (std::vector<size_t> const& previous_class_labels, std::vector<size_t> const& previous_cluster_labels, std::vector<size_t> const& current_class_labels, std::vector<size_t> const& current_cluster_labels)
    {
        // check if the shape of some partitions didn't change
        // for this do parallel DFS on the cluster labels from previous and current iteration
        // partitions do not change iff their DFS traces are the exactly the same

        auto const number_of_components = *std::max_element(current_cluster_labels.begin(), current_cluster_labels.end()) + 1;

        std::vector<char> changed(number_of_components);

        std::stack<size_t> S;
        std::vector<char> visited(problem.numberOfVertices());
        
        std::vector<size_t> vertices_of_cluster;
        for (size_t i = 0; i < problem.numberOfVertices(); ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto const label_new = current_cluster_labels[i];
                auto const label_old = previous_cluster_labels[i];
                
                while (!S.empty())
                {
                    auto const v = S.top();
                    S.pop();

                    vertices_of_cluster.push_back(v);

                    for (auto w = problem.originalGraph().verticesFromVertexBegin(v); w != problem.originalGraph().verticesFromVertexEnd(v); ++w)
                    {
                        if (
                            previous_cluster_labels[*w] == label_old && current_cluster_labels[*w] != label_new
                            ||
                            previous_cluster_labels[*w] != label_old && current_cluster_labels[*w] == label_new
                            )
                            changed[label_new] = 1;

                        if (visited[*w])
                            continue;

                        if (current_cluster_labels[*w] == label_new)
                        {
                            S.push(*w);
                            visited[*w] = 1;

                            if (previous_cluster_labels[*w] != label_old)
                                changed[label_new] = 1;
                        }
                    }
                }

                if (problem.liftedGraph().numberOfEdges() > problem.originalGraph().numberOfEdges())
                    for (auto v : vertices_of_cluster)
                    {
                        if (current_class_labels[v] != previous_class_labels[v])
                        {
                            changed[current_cluster_labels[v]] = 1;
                            break;
                        }

                        for (auto w = problem.liftedGraph().verticesFromVertexBegin(v); w != problem.liftedGraph().verticesFromVertexEnd(v); ++w)
                            if (current_class_labels[*w] != previous_class_labels[*w])
                            {
                                changed[current_cluster_labels[v]] = 1;
                                break;
                            }

                        if (changed[current_cluster_labels[v]])
                            break;
                    }

                vertices_of_cluster.clear();
            }

        return changed;
    };

    auto project_onto_feasable_set = [&problem] (std::vector<size_t> const& cluster_labels)
    {
        std::stack<size_t> S;
        std::vector<char> visited(problem.numberOfVertices());

        std::vector<size_t> feasable_cluster_labels(cluster_labels.size());

        // do connected component labeling on the original graph
        for (size_t i = 0, new_label = 0; i < problem.numberOfVertices(); ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto const label = cluster_labels[i];

                while (!S.empty())
                {
                    auto const v = S.top();
                    S.pop();

                    feasable_cluster_labels[v] = new_label;

                    for (auto it = problem.originalGraph().adjacenciesFromVertexBegin(v); it != problem.originalGraph().adjacenciesFromVertexEnd(v); ++it)
                        if (cluster_labels[it->vertex()] == label && !visited[it->vertex()])
                        {
                            S.push(it->vertex());
                            visited[it->vertex()] = 1;
                        }
                }

                ++new_label;
            }

        return feasable_cluster_labels;
    };


    for (size_t i = 0; i < problem.numberOfVertices(); ++i)
    {
        buffer.class_labels[i] = input[i].classIndex;
        buffer.cluster_labels[i] = input[i].clusterIndex;
    }

    buffer.cluster_labels = project_onto_feasable_set(buffer.cluster_labels);

    auto partitions = form_partitions(buffer.cluster_labels);
    buffer.max_not_used_label = partitions.size();

    double current_obj_value = compute_obj_value(buffer.class_labels, buffer.cluster_labels);

    std::cout << "starting objective: " << std::fixed << current_obj_value << std::endl;

    // we might need to rollback to a previous solution, so keep it
    auto last_good_class_labels = buffer.class_labels;
    auto last_good_cluster_labels = buffer.cluster_labels;

    auto class_labels_for_changed_thing = buffer.class_labels;
    auto cluster_labels_for_changed_thing = buffer.cluster_labels;

    // interatively update bipartition in order to minimize the total cost of the multicut
    for (size_t k = 0; ; ++k)
    {
        // 1 if i-th partitioned changed since last iteration, 0 otherwise
        std::vector<char> changed(partitions.size(), 1);

        if (k > 0)
            changed = mark_partitions_that_changed(class_labels_for_changed_thing, cluster_labels_for_changed_thing, buffer.class_labels, buffer.cluster_labels);
        
        // build partions' adjacency graph as given by original graph connectivity
        std::vector<std::set<size_t>> edges(partitions.size());
        for (size_t e = 0; e < problem.originalGraph().numberOfEdges(); ++e)
        {
            auto const v0 = last_good_cluster_labels[problem.originalGraph().vertexOfEdge(e, 0)];
            auto const v1 = last_good_cluster_labels[problem.originalGraph().vertexOfEdge(e, 1)];

            if (v0 != v1)
                edges[std::min(v0, v1)].insert(std::max(v0, v1));
        }

        // update boundary between all neighbors in the above computed partition graph
        // skip pairs if both partitions didn't change since last time
        for (size_t i = 0; i < partitions.size(); ++i)
            if (!partitions[i].empty())
                for (auto j : edges[i])
                    if (!partitions[j].empty() && (changed[j] || changed[i]))
                    {
                        update_bipartition(partitions[i], partitions[j]);

                        if (partitions[i].size() == 0)
                            break;
                    }
        
        
        buffer.cluster_labels = project_onto_feasable_set(buffer.cluster_labels);

        partitions = form_partitions(buffer.cluster_labels);
        buffer.max_not_used_label = partitions.size();

        auto pair_updates_decrease = current_obj_value - compute_obj_value(buffer.class_labels, buffer.cluster_labels);
        current_obj_value -= pair_updates_decrease;


        if (k > 0)
            changed = mark_partitions_that_changed(class_labels_for_changed_thing, cluster_labels_for_changed_thing, buffer.class_labels, buffer.cluster_labels);
        else
            changed.resize(partitions.size(), 1);

        // try to intoduce new partitions
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

        buffer.cluster_labels = project_onto_feasable_set(buffer.cluster_labels);

        auto new_sets_decrease = current_obj_value - compute_obj_value(buffer.class_labels, buffer.cluster_labels);
        current_obj_value -= new_sets_decrease;

        auto update_labels_decrease = .0;
        if (problem.numberOfClasses() > 1)
        {
            std::map<size_t, size_t> indices;
            for (size_t v = 0; v < problem.numberOfVertices(); ++v)
                indices[buffer.cluster_labels[v]] = 0;

            size_t cnt = 0;
            for (auto& p : indices)
                p.second = cnt++;

            Solution tmp(problem.numberOfVertices());
            for (size_t v = 0; v < problem.numberOfVertices(); ++v)
            {
                tmp[v].clusterIndex = indices[buffer.cluster_labels[v]];
                tmp[v].classIndex = buffer.class_labels[v];
            }

            update_labels_decrease = detail::update_class_labels(problem, tmp, tmp);

            current_obj_value -= update_labels_decrease;

            for (size_t v = 0; v < problem.numberOfVertices(); ++v)
            {
                buffer.cluster_labels[v] = tmp[v].clusterIndex;
                buffer.class_labels[v] = tmp[v].classIndex;
            }
        }        
        
        // if the new true energy is higher, than the current one, terminate
        if (pair_updates_decrease + new_sets_decrease + update_labels_decrease < 1e-6)
            break;

        partitions = form_partitions(buffer.cluster_labels);
        buffer.max_not_used_label = partitions.size();

        // in order to avoid possible numerical instability check whether the solution has actually changed since last time
        bool didnt_change = true;
        for (size_t i = 0; i < problem.originalGraph().numberOfEdges(); ++i)
        {
            auto v0 = problem.originalGraph().vertexOfEdge(i, 0);
            auto v1 = problem.originalGraph().vertexOfEdge(i, 1);

            if ((buffer.cluster_labels[v0] != buffer.cluster_labels[v1]) != (last_good_cluster_labels[v0] != last_good_cluster_labels[v1]))
            {
                didnt_change = false;
                break;
            }
        }

        for (size_t v = 0; v < problem.originalGraph().numberOfVertices(); ++v)
            if (last_good_class_labels[v] != buffer.class_labels[v])
            {
                didnt_change = false;
                break;
            }

        if (didnt_change)
            break;

        visitor(buffer.cluster_labels, buffer.class_labels);

        std::cout << "....pair updates: " << pair_updates_decrease << std::endl;
        std::cout << "....new sets: " << new_sets_decrease << std::endl;

        if (problem.numberOfClasses() > 1)
            std::cout << "....updating class labels: " << update_labels_decrease << std::endl;
        
        std::cout << "..new objective: " << current_obj_value << std::endl;
        
        class_labels_for_changed_thing = last_good_class_labels;
        cluster_labels_for_changed_thing = last_good_cluster_labels;

        // if everything is fine with the new solution, remember it
        last_good_class_labels = buffer.class_labels;
        last_good_cluster_labels = buffer.cluster_labels;
    }

    // the following simply makes sure that the cluster labels are in [0, number_of_clusters)
    std::map<size_t, size_t> indices;
    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
        indices[last_good_cluster_labels[v]] = 0;

    size_t cnt = 0;
    for (auto& p : indices)
        p.second = cnt++;

    Solution output(problem.numberOfVertices());

    for (size_t v = 0; v < problem.numberOfVertices(); ++v)
    {
        output[v].clusterIndex = indices[last_good_cluster_labels[v]];
        output[v].classIndex = last_good_class_labels[v];
    }

    visitor(output);

    return output;
}

}

#endif
