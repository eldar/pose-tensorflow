# distutils: language = c++

cimport cython
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp cimport bool

np.import_array()

cdef extern from "solve_nl_lmp.hxx":
  void solve_nl_lmp_cpp(double* unValues, int, int,
                        unsigned short* pwIndices, int, int,
                        double *pwValues, int, int,
                        bool is_sparse_graph, bool solver_type, bool do_suppression, bool do_logit_transform,
                        unsigned long long *result)


@cython.boundscheck(False)
@cython.wraparound(False)
def solve_nl_lmp(np.ndarray[np.float64_t, ndim=2, mode="c"] unary_array,
                 np.ndarray[np.uint16_t, ndim=2, mode="c"] pwidx_array,
                 np.ndarray[np.float64_t, ndim=2, mode="c"] pw_array,
                 is_sparse_graph, solver_type, do_suppression, logit_in_solver):


  cdef np.ndarray[np.uint64_t, ndim=2, mode="c"] result = np.zeros([unary_array.shape[0], 2], dtype=np.uint64)

  solve_nl_lmp_cpp(&unary_array[0, 0], unary_array.shape[0], unary_array.shape[1],
                   &pwidx_array[0, 0], pwidx_array.shape[0], pwidx_array.shape[1],
                   &pw_array[0, 0], pw_array.shape[0], pw_array.shape[1],
                   is_sparse_graph, solver_type, do_suppression, logit_in_solver,
                   &result[0, 0])

  return result
