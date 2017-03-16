/// \file
/// \brief The gsl_matrix_wrp class header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_GSL_GSL_MATRIX_WRP_H
#define _CATH_TOOLS_SOURCE_COMMON_GSL_GSL_MATRIX_WRP_H

#include "exception/runtime_error_exception.hpp"

#include <gsl/gsl_matrix.h>

namespace cath {
	namespace geom {
		namespace detail {

			/// \brief RAII wrapper for a GSL gsl_matrix pointer
			class gsl_matrix_wrp final {
			private:
				/// \brief Pointer to a GSL gsl_matrix
				gsl_matrix * ptr;

			public:
				/// \brief Delete the default ctor
				gsl_matrix_wrp() = delete;

				/// \brief Ctor that wraps call to gsl_matrix_alloc()
				gsl_matrix_wrp(const size_t &arg_n1, ///< The first  dimension of the matrix
				               const size_t &arg_n2  ///< The second dimension of the matrix
				               ) : ptr{ gsl_matrix_alloc( arg_n1, arg_n2 ) } {
					if ( ptr == nullptr ) {
						BOOST_THROW_EXCEPTION(common::runtime_error_exception("Was unable to allocate gsl_matrix"));
					}
				}

				/// \brief Dtor to wrap gsl_matrix_free()
				~gsl_matrix_wrp() noexcept {
					try {
						gsl_matrix_free( ptr );
					}
					catch (...) {
					}
				}

				/// \brief Const-overload of getter for a reference to the gsl_matrix
				const gsl_matrix & get_ref() const {
					return *ptr;
				}

				/// \brief Non-const-overload of getter for a reference to the gsl_matrix
				gsl_matrix & get_ref() {
					return *ptr;
				}

				/// \brief Const-overload of getter for a pointer to the gsl_matrix
				const gsl_matrix * get_ptr() const {
					return ptr;
				}

				/// \brief Non-const-overload of getter for a pointer to the gsl_matrix
				gsl_matrix * get_ptr() {
					return ptr;
				}
			};

		} // namespace detail
	} // namespace geom
} // namespace cath

#endif
