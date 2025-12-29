#include "vec3.h"
#include "pack.h"

#include "independent.h"
#include "jacobian.h"
#include "schemes/include.h"
#include "boundary_manager.h"
#include "utils/Timer.h"

#include "projected_saver_2D.hpp"
#include "parameters.h"

#include "FDTD/fdtd_solver_1d.h"
#include "calc_current_in_x_and_y.h"

#include "normalization.h"

#include "none.h"
#include "config.h"

#ifdef USE_SUPERCOMPUTER
#include "supercomputer_instruments/advection_equation.h"
#include "supercomputer_instruments/axis.h"
#include "supercomputer_instruments/n_d_tensor_with_ghost_cell.h"
#include "supercomputer_instruments/n_d_tensor.h"
#else
#include "advection_equation.h"
#include "axis.h"
#include "n_d_tensor_with_ghost_cell.h"
#include "n_d_tensor.h"
#endif
