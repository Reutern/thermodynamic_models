/* siman/gsl_siman.h
*
* Copyright (C) 1996, 1997, 1998, 1999, 2000 Mark Galassi
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/


#ifndef __GSL_SIMAN_H__
#define __GSL_SIMAN_H__
#include <stdlib.h>
#include </usr/local/include/gsl/gsl_rng.h>
#include </usr/local/include/gsl/gsl_math.h>
#include </usr/local/include/gsl/gsl_multimin.h>
#include <iostream>

/* types for the function pointers passed to gsl_siman_solve */

typedef double (*gsl_siman_Efunc_t) (const gsl_vector* xp, void* par);
typedef void (*gsl_siman_step_t) (const gsl_rng *r, gsl_vector* xp, double step_size);
typedef double (*gsl_siman_metric_t) (gsl_vector* xp, gsl_vector* yp);
typedef void (*gsl_siman_print_t) (gsl_vector* xp);
typedef void (*gsl_siman_copy_t) (gsl_vector *source, gsl_vector *dest);
typedef void * (*gsl_siman_copy_construct_t) (gsl_vector* xp);
typedef void (*gsl_siman_destroy_t) (gsl_vector* xp);

/* this structure contains all the information needed to structure the
search, beyond the energy function, the step function and the
initial guess. */

typedef struct {
  int n_tries; /* how many points to try for each step */
  int iters_fixed_T; /* how many iterations at each temperature? */
  double step_size; /* max step size in the random walk */
  /* the following parameters are for the Boltzmann distribution */
  double k, t_initial, mu_t, t_min;
} gsl_siman_params_t;

/* prototype for the workhorse function */

void gsl_siman_solve(const gsl_rng * r,
                     gsl_vector *x0_p, void* predictor , gsl_siman_Efunc_t Ef,
                     gsl_siman_step_t take_step,
                     gsl_siman_print_t print_position,
                     gsl_siman_copy_t copyfunc,
                     gsl_siman_copy_construct_t copy_constructor,
                     gsl_siman_destroy_t destructor,
                     size_t element_size,
                     gsl_siman_params_t params);

void
gsl_siman_solve_many (const gsl_rng * r, gsl_vector *x0_p, void* predictor, gsl_siman_Efunc_t Ef,
                      gsl_siman_step_t take_step,
                      gsl_siman_metric_t distance,
                      gsl_siman_print_t print_position,
                      size_t element_size,
                      gsl_siman_params_t params);

#endif /* __GSL_SIMAN_H__ */
