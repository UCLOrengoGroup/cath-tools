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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GSL_VECTOR_WRP_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GSL_VECTOR_WRP_HPP

#include "cath/common/exception/runtime_error_exception.hpp"

#include <gsl/gsl_vector.h>

namespace cath::geom::detail {

	/// \brief RAII wrapper for a GSL gsl_vector pointer
	class gsl_vector_wrp final {
	private:
		/// \brief Pointer to a GSL gsl_vector
		gsl_vector * ptr;

	public:
		/// \brief Delete the default ctor
		gsl_vector_wrp() = delete;

		/// \brief Ctor that wraps call to gsl_vector_alloc()
		explicit gsl_vector_wrp(const size_t &prm_n ///< The dimension of the vector
		                        ) : ptr{ gsl_vector_alloc( prm_n ) } {

			if ( ptr == nullptr ) {
				BOOST_THROW_EXCEPTION(common::runtime_error_exception("Was unable to allocate gsl_vector"));
			}
		}

		/// \brief Dtor to wrap gsl_vector_free()
		~gsl_vector_wrp() noexcept {
			try {
				if ( ptr != nullptr ) {
					gsl_vector_free( ptr );
				}
			}
			catch (...) {
			}
		}

		gsl_vector_wrp(const gsl_vector_wrp &) = delete; ///< Make move-only

		/// \brief Move ctor that sets RHS's ptr to null after the move
		gsl_vector_wrp(gsl_vector_wrp &&prm_rhs ///< The gsl_vector_wrp to move into this
		               ) : ptr { std::move( prm_rhs.ptr ) } {
			prm_rhs.ptr = nullptr;
		}

		gsl_vector_wrp & operator=(const gsl_vector_wrp &) = delete; ///< Make move-only

		/// \brief Move assignment operator that sets RHS's ptr to null after the move
		gsl_vector_wrp & operator=(gsl_vector_wrp &&prm_rhs ///< The gsl_vector_wrp to move into this
		                           ) {
			if ( this != &prm_rhs ) {
				ptr = std::move( prm_rhs.ptr );
				prm_rhs.ptr = nullptr;
			}
			return *this;
		}

		/// \brief Const-overload of getter for a reference to the gsl_vector
		[[nodiscard]] const gsl_vector &get_ref() const {
			return *ptr;
		}

		/// \brief Non-const-overload of getter for a reference to the gsl_vector
		gsl_vector & get_ref() {
			return *ptr;
		}

		/// \brief Const-overload of getter for a pointer to the gsl_vector
		[[nodiscard]] const gsl_vector *get_ptr() const {
			return ptr;
		}

		/// \brief Non-const-overload of getter for a pointer to the gsl_vector
		gsl_vector * get_ptr() {
			return ptr;
		}
	};

	/// \brief Generate a string describing the specified gsl_vector_wrp
	///
	/// \relates to_string
	inline std::string to_string(const gsl_vector_wrp &prm_vector_wrap ///< The gsl_vector_wrp to describe
	                             ) {
		using ::std::to_string;
		return
			  "gsl_vector_wrp["
			+ to_string( gsl_vector_get( prm_vector_wrap.get_ptr(), 0 ) )
			+ ", "
			+ to_string( gsl_vector_get( prm_vector_wrap.get_ptr(), 1 ) )
			+ ", "
			+ to_string( gsl_vector_get( prm_vector_wrap.get_ptr(), 2 ) )
			+ "]";
	}


} // namespace cath::geom::detail

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_GSL_GSL_VECTOR_WRP_HPP
