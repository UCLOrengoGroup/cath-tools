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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GSL_MATRIX_WRP_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GSL_MATRIX_WRP_HPP

#include "cath/common/exception/runtime_error_exception.hpp"

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
				gsl_matrix_wrp(const size_t &prm_n1, ///< The first  dimension of the matrix
				               const size_t &prm_n2  ///< The second dimension of the matrix
				               ) : ptr{ gsl_matrix_alloc( prm_n1, prm_n2 ) } {
					if ( ptr == nullptr ) {
						BOOST_THROW_EXCEPTION(common::runtime_error_exception("Was unable to allocate gsl_matrix"));
					}
				}

				/// \brief Dtor to wrap gsl_matrix_free()
				~gsl_matrix_wrp() noexcept {
					try {
						if ( ptr != nullptr ) {
							gsl_matrix_free( ptr );
						}
					}
					catch (...) {
					}
				}

				gsl_matrix_wrp(const gsl_matrix_wrp &) = delete; ///< Make move-only

				/// \brief Move ctor that sets RHS's ptr to null after the move
				gsl_matrix_wrp(gsl_matrix_wrp &&prm_rhs ///< The gsl_matrix_wrp to move into this
				               ) : ptr { std::move( prm_rhs.ptr ) } {
					prm_rhs.ptr = nullptr;
				}

				gsl_matrix_wrp & operator=(const gsl_matrix_wrp &) = delete; ///< Make move-only

				/// \brief Move assignment operator that sets RHS's ptr to null after the move
				gsl_matrix_wrp & operator=(gsl_matrix_wrp &&prm_rhs ///< The gsl_matrix_wrp to move into this
				                           ) {
					if ( this != &prm_rhs ) {
						ptr = std::move( prm_rhs.ptr );
						prm_rhs.ptr = nullptr;
					}
					return *this;
				}

				/// \brief Const-overload of getter for a reference to the gsl_matrix
				[[nodiscard]] const gsl_matrix &get_ref() const {
					return *ptr;
				}

				/// \brief Non-const-overload of getter for a reference to the gsl_matrix
				gsl_matrix & get_ref() {
					return *ptr;
				}

				/// \brief Const-overload of getter for a pointer to the gsl_matrix
				[[nodiscard]] const gsl_matrix *get_ptr() const {
					return ptr;
				}

				/// \brief Non-const-overload of getter for a pointer to the gsl_matrix
				gsl_matrix * get_ptr() {
					return ptr;
				}
			};

			/// \brief Generate a string describing the specified gsl_matrix_wrp
			///
			/// \relates to_string
			inline std::string to_string(const gsl_matrix_wrp &prm_matrix_wrap ///< The gsl_matrix_wrp to describe
			                             ) {
				using ::std::to_string;
				return "gsl_matrix_wrp["
					       + to_string( gsl_matrix_get( prm_matrix_wrap.get_ptr(), 0, 0 ) )
					+ ", " + to_string( gsl_matrix_get( prm_matrix_wrap.get_ptr(), 0, 1 ) )
					+ ", " + to_string( gsl_matrix_get( prm_matrix_wrap.get_ptr(), 0, 2 ) )

					+ ", " + to_string( gsl_matrix_get( prm_matrix_wrap.get_ptr(), 1, 0 ) )
					+ ", " + to_string( gsl_matrix_get( prm_matrix_wrap.get_ptr(), 1, 1 ) )
					+ ", " + to_string( gsl_matrix_get( prm_matrix_wrap.get_ptr(), 1, 2 ) )

					+ ", " + to_string( gsl_matrix_get( prm_matrix_wrap.get_ptr(), 2, 0 ) )
					+ ", " + to_string( gsl_matrix_get( prm_matrix_wrap.get_ptr(), 2, 1 ) )
					+ ", " + to_string( gsl_matrix_get( prm_matrix_wrap.get_ptr(), 2, 2 ) )
					+ "]";
			}

			/// \brief Increment the specified cell of the specified gsl_matrix_wrp by the specified amount
			///
			/// \relates gsl_matrix_wrp
			inline void gsl_matrix_wrp_increment(gsl_matrix_wrp &prm_matrix_wrap, ///< The gsl_matrix_wrp to modify
			                                     const size_t   &prm_row_index,   ///< The row_index of the cell to modify
			                                     const size_t   &prm_col_index,   ///< The column index of the cell to modify
			                                     const double   &prm_addend       ///< The amount to add to the cell
			                                     ) {
				gsl_matrix_set(
					prm_matrix_wrap.get_ptr(),
					prm_row_index,
					prm_col_index,
					gsl_matrix_get(
						prm_matrix_wrap.get_ptr(),
						prm_row_index,
						prm_col_index
					) + prm_addend
				);
			}

		} // namespace detail
	} // namespace geom
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GSL_MATRIX_WRP_HPP
