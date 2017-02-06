/// \file
/// \brief The eigen class header

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
///
/// This program is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _CATH_TOOLS_SOURCE_EIGEN_EIGEN_H
#define _CATH_TOOLS_SOURCE_EIGEN_EIGEN_H

// #include "eigen/eigen_type_aliases.hpp"
// #include "eigen/region/region.hpp"

#include <cstddef>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

namespace cath {
	namespace chop {

		int my_function() {
			double data[] = { 1.0  , 1/2.0, 1/3.0, 1/4.0,
			                  1/2.0, 1/3.0, 1/4.0, 1/5.0,
			                  1/3.0, 1/4.0, 1/5.0, 1/6.0,
			                  1/4.0, 1/5.0, 1/6.0, 1/7.0 };

			gsl_matrix_view m  = gsl_matrix_view_array (data, 4, 4);

			gsl_vector *eval = gsl_vector_alloc (4);
			gsl_matrix *evec = gsl_matrix_alloc (4, 4);

			gsl_eigen_symmv_workspace * w =  gsl_eigen_symmv_alloc (4);

			gsl_eigen_symmv (&m.matrix, eval, evec, w);

			gsl_eigen_symmv_free (w);

			gsl_eigen_symmv_sort (eval, evec,  GSL_EIGEN_SORT_ABS_ASC);

			for (size_t i = 0; i < 4; ++i) {
				double          eval_i  = gsl_vector_get   ( eval, i );
				gsl_vector_view evec_i  = gsl_matrix_column( evec, i );

				printf( "eigenvalue = %g\n", eval_i );
				printf( "eigenvector = \n" );
				gsl_vector_fprintf( stdout,  &evec_i.vector, "%g" );
			}

			gsl_vector_free( eval );
			gsl_matrix_free( evec );

			return 0;
		}

	} // namespace chop
} // namespace cath

#endif
