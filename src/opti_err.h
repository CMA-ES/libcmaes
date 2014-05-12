/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 Inria
 * Author: Emmanuel Benazera <emmanuel.benazera@lri.fr>
 *
 * This file is part of libcmaes.
 *
 * libcmaes is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libcmaes is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libcmaes.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OPTI_ERR_H
#define OPTI_ERR_H

/* errors */
enum {
  /* reaches convergence. */
  OPTI_SUCCESS = 0,
  OPTI_STOP,
  /* intial variables already minimize the objective function. */
  OPTI_ALREADY_MINIMIZED,
  /* unknown error */
  OPTI_ERR_UNKNOWN = -1024,
  /* insufficient memory. */
  OPTI_ERR_OUTOFMEMORY=-1025,
  /* invalid number of variables specified. */
  OPTI_ERR_INVALID_N=-1026,
  /* the algorithm has reached a termination criteria without reaching the objective. */
  OPTI_ERR_TERMINATION=-1
};

#endif
