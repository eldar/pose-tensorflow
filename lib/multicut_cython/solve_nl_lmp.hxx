#define HAVE_CPP11_INITIALIZER_LISTS
#include "andres/marray.hxx"

#include "nl-lmp.hxx"

#include <algorithm>

#include <vector>

#include <iostream>

using namespace std;

template <typename T> 
andres::View<T> get_view2d(T *gridData, int grid_height, int grid_width)
{
  std::array<size_t, 2> shape;
  shape[0] = grid_height;
  shape[1] = grid_width;
  andres::View<T> view(shape.begin(), shape.end(), gridData, andres::FirstMajorOrder); 
  return view;
}


void solve_nl_lmp_cpp(double *unValues, int un_H, int un_W,
                      uint16_t* pwIndices, int pwi_H, int pwi_W,
                      double* pwValues, int pwv_H, int pwv_W,
                      bool is_sparse_graph, bool solver_type, bool do_suppression, bool do_logit_transform,
                      uint64_t *result)
{
  andres::View<uint16_t> pwind_view = get_view2d<uint16_t>(pwIndices, pwi_H, pwi_W);
  andres::View<double> un_view = get_view2d<double>(unValues, un_H, un_W);
  andres::View<double> pw_view = get_view2d<double>(pwValues, pwv_H, pwv_W);

  const bool is_node_labeling_problem = (pw_view.shape(1) == 4);

  andres::Marray<size_t> unLab;

  SolvingMethod method = SolvingMethod::Joint;

  if (is_node_labeling_problem){
    // SolvingMethod method;
    // if(solver_type) {
    //   method = SolvingMethod::Joint;
    //   std::cout << "Method: Joint" << std::endl;
    // }
    // else {
    //   method = SolvingMethod::Simple;
    //   std::cout << "Method: Simple" << std::endl;
    // }

    assert(false);
  }
  else {
    if (is_sparse_graph) {
      unLab = solve_nl_lmp_sparse_graph(un_view, pwind_view, pw_view, method, do_logit_transform, do_suppression);
    }
    else {
      //unLab = solve_lmp_complete_graph(un_view, pwind_view, pw_view, do_suppression);

      unLab = solve_nl_lmp_complete_graph(un_view, pwind_view, pw_view, method, do_logit_transform, do_suppression);  
    }
  }
  assert(unLab.size() > 0);

  andres::View<uint64_t> result_view = get_view2d<uint64_t>(result, un_H, 2);
  for(int k = 0; k < un_H; ++k)
  {
    result_view(k, 0) = unLab(k, 0);
    result_view(k, 1) = unLab(k, 1);
  }
}
