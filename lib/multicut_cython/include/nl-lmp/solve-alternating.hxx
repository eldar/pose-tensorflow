#pragma once
#ifndef NL_LMP_SOLVE_ALTERNATING_HXX
#define NL_LMP_SOLVE_ALTERNATING_HXX

#include "detail/update-class-labels.hxx"
#include "detail/call-multicut-solver.hxx"



namespace nl_lmp
{

template<typename GRAPH>
inline
Solution solve_alternating(Problem<GRAPH> const& problem, Solution const& input)
{
    struct Visitor
    {
        void operator()(Solution const& solution) const
        {}
    } visitor;

    return solve_alternating(problem, input, visitor);
}

template<typename GRAPH, typename VISITOR>
inline
Solution solve_alternating(Problem<GRAPH> const& problem, Solution const& input, VISITOR& visitor)
{
    double obj_value = detail::compute_obj_value(problem, input);

    Solution output(input);

    std::cout << "starting objective value: " << std::fixed << obj_value << std::endl;

    double gain_update_labels = .0;
    double gain_update_multicut = .0;
    do
    {
        gain_update_multicut = detail::update_multicut(problem, output, output);

        gain_update_labels = detail::update_class_labels(problem, output, output);

        visitor(output);

        if (gain_update_labels + gain_update_multicut > 1e-6)
        {
            obj_value -= gain_update_labels + gain_update_multicut;

            std::cout << "..gain from updating multicut: " << gain_update_multicut << std::endl;
            std::cout << "..gain from updating labels: " << gain_update_labels << std::endl;
            std::cout << "..new objective value: " << obj_value << std::endl;
        }
    } while (gain_update_labels + gain_update_multicut > 1e-6);

    return output;
}

}

#endif
