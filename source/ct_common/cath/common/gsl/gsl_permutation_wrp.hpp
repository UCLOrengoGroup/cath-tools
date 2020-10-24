/// \file
/// \brief The gsl_permutation_wrp class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GSL_PERMUTATION_WRP_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GSL_PERMUTATION_WRP_HPP

#include "cath/common/exception/runtime_error_exception.hpp"

#include <gsl/gsl_permutation.h>

namespace cath {
	namespace geom {
		namespace detail {

			/// \brief RAII wrapper for a GSL gsl_permutation pointer
			class gsl_permutation_wrp final {
			private:
				/// \brief Pointer to a GSL gsl_permutation
				gsl_permutation * ptr;

			public:
				/// \brief Delete the default ctor
				gsl_permutation_wrp() = delete;

				/// \brief Ctor that wraps call to gsl_permutation_calloc() (which allocates *and* initialises)
				gsl_permutation_wrp(const size_t &prm_n ///< The dimension of the permutation
				                    ) : ptr{ gsl_permutation_calloc( prm_n ) } {
					if ( ptr == nullptr ) {
						BOOST_THROW_EXCEPTION(common::runtime_error_exception("Was unable to allocate+initialise gsl_permutation"));
					}
				}

				/// \brief Dtor to wrap gsl_permutation_free()
				~gsl_permutation_wrp() noexcept {
					try {
						if ( ptr != nullptr ) {
							gsl_permutation_free( ptr );
						}
					}
					catch (...) {
					}
				}

				gsl_permutation_wrp(const gsl_permutation_wrp &) = delete; ///< Make move-only

				/// \brief Move ctor that sets RHS's ptr to null after the move
				gsl_permutation_wrp(gsl_permutation_wrp &&prm_rhs ///< The gsl_permutation_wrp to move into this
				                    ) : ptr { std::move( prm_rhs.ptr ) } {
					prm_rhs.ptr = nullptr;
				}

				gsl_permutation_wrp & operator=(const gsl_permutation_wrp &) = delete; ///< Make move-only

				/// \brief Move assignment operator that sets RHS's ptr to null after the move
				gsl_permutation_wrp & operator=(gsl_permutation_wrp &&prm_rhs ///< The gsl_permutation_wrp to move into this
				                           ) {
					if ( this != &prm_rhs ) {
						ptr = std::move( prm_rhs.ptr );
						prm_rhs.ptr = nullptr;
					}
					return *this;
				}

				/// \brief Const-overload of getter for a reference to the gsl_permutation
				const gsl_permutation & get_ref() const {
					return *ptr;
				}

				/// \brief Non-const-overload of getter for a reference to the gsl_permutation
				gsl_permutation & get_ref() {
					return *ptr;
				}

				/// \brief Const-overload of getter for a pointer to the gsl_permutation
				const gsl_permutation * get_ptr() const {
					return ptr;
				}

				/// \brief Non-const-overload of getter for a pointer to the gsl_permutation
				gsl_permutation * get_ptr() {
					return ptr;
				}
			};

		} // namespace detail
	} // namespace geom
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GSL_PERMUTATION_WRP_HPP
