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

#ifndef LLOGGING_H
#define LLOGGING_H

#include "libcmaes_config.h"

#ifdef HAVE_GLOG // HAVE_LIB_GLOG
#include <glog/logging.h>
#else
#include <iostream>

namespace libcmaes
{
  static std::string INFO="INFO";
  static std::string WARNING="WARNING";
  static std::string ERROR="ERROR";
  static std::string FATAL="FATAL";

  static std::ostream nullstream(0);

inline std::ostream& LOG(const std::string &severity,std::ostream &out=std::cout)
{
  out << severity << " - ";
  return out;
}

inline std::ostream& LOG_IF(const std::string &severity,const bool &condition,std::ostream &out=std::cout)
{
  if (condition)
    return LOG(severity,out);
  else return nullstream;
}
}
#endif
#endif
