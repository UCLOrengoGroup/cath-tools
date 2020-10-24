/// \file
/// \brief The get_determinant class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GET_DETERMINANT_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GET_DETERMINANT_HPP

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/gsl/gsl_matrix_wrp.hpp"
#include "cath/common/gsl/gsl_permutation_wrp.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

namespace cath {
	namespace geom {
		namespace detail {

			/// \brief Get the determinant of the specified gsl_matrix_wrp
			inline double get_determinant(const gsl_matrix_wrp &prm_matrix ///< The gsl_matrix_wrp for which the determinant should be calculated
				                          ) {
				int          sign   = 1;
				const auto   size   = prm_matrix.get_ref().size1;
				if ( size != prm_matrix.get_ref().size2 ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot calculate determinant of a matrix that isn't square"));
				}
				gsl_matrix_wrp      temp_matrix {  size, size };
				gsl_permutation_wrp temp_permutn{ size        };

				gsl_matrix_memcpy   ( temp_matrix.get_ptr(), prm_matrix.get_ptr()          );
				gsl_linalg_LU_decomp( temp_matrix.get_ptr(), temp_permutn.get_ptr(), &sign );

				return gsl_linalg_LU_det( temp_matrix.get_ptr(), sign );
			}

		} // namespace detail
	} // namespace geom
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GET_DETERMINANT_HPP
