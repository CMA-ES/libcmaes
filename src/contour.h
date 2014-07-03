/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 INRIA
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

#ifndef CONTOUR_H
#define CONTOUR_H

#include <vector>

namespace libcmaes
{

  /**
   * \brief function contour as a set of points and values.
   */
  class contour
  {
  public:
    contour() {}
    ~contour() {}

    /**
     * \brief add a contour point.
     * @param x parameter value
     * @param y parameter value
     */
    void add_point(const double &x, const double &y)
    {
      _points.push_back(std::pair<double,double>(x,y));
    }

    std::ostream& print(std::ostream &out) const
      {
	out << "contour points: [";
	for (size_t i=0;i<_points.size();i++)
	  {
	    out << "(" << _points.at(i).first << "," << _points.at(i).second << ") ";
	  }
	out << "]\n";
	return out;
      }
    
    std::vector<std::pair<double,double>> _points;
  };

  std::ostream& operator<<(std::ostream &out, const contour &c)
    {
      c.print(out);
      return out;
    }
  
}

#endif
