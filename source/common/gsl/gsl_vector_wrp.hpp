/// \file
/// \brief The gsl_vector_wrp class header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_GSL_GSL_VECTOR_WRP_H
#define _CATH_TOOLS_SOURCE_COMMON_GSL_GSL_VECTOR_WRP_H

#include "exception/runtime_error_exception.hpp"

#include <gsl/gsl_vector.h>

namespace cath {
	namespace geom {
		namespace detail {

			/// \brief RAII wrapper for a GSL gsl_vector pointer
			class gsl_vector_wrp final {
			private:
				/// \brief Pointer to a GSL gsl_vector
				gsl_vector * ptr;

			public:
				/// \brief Delete the default ctor
				gsl_vector_wrp() = delete;

				/// \brief Ctor that wraps call to gsl_vector_alloc()
				gsl_vector_wrp(const size_t &arg_n ///< The dimension of the vector
				               ) : ptr{ gsl_vector_alloc( arg_n ) } {

					if ( ptr == nullptr ) {
						BOOST_THROW_EXCEPTION(common::runtime_error_exception("Was unable to allocate gsl_vector"));
					}
				}

				/// \brief Dtor to wrap gsl_vector_free()
				~gsl_vector_wrp() noexcept {
					try {
						gsl_vector_free( ptr );
					}
					catch (...) {
					}
				}

				/// \brief Const-overload of getter for a reference to the gsl_vector
				const gsl_vector & get_ref() const {
					return *ptr;
				}

				/// \brief Non-const-overload of getter for a reference to the gsl_vector
				gsl_vector & get_ref() {
					return *ptr;
				}

				/// \brief Const-overload of getter for a pointer to the gsl_vector
				const gsl_vector * get_ptr() const {
					return ptr;
				}

				/// \brief Non-const-overload of getter for a pointer to the gsl_vector
				gsl_vector * get_ptr() {
					return ptr;
				}
			};

		} // namespace detail
	} // namespace geom
} // namespace cath

#endif
