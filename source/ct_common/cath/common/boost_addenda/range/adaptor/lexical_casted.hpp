/// \file
/// \brief The lexical_casted class header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_LEXICAL_CASTED_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_LEXICAL_CASTED_HPP

#include "cath/common/boost_addenda/range/adaptor/range/lexical_casted_range.hpp"

namespace cath::common {

	/// \brief TODOCUMENT
	template <typename T>
	class lexical_casted final {};

	/// \brief Non-const range overload of operator| for lexical_casted range adaptor
	template <typename T, typename ForwardRng>
	inline detail::lexical_casted_range<T, ForwardRng> operator|(ForwardRng              &prm_range, ///< The range to which the lexical_casted adaptor should be applied
	                                                             const lexical_casted<T> &           ///< A dummy lexical_casted parameter for determining which adaptor should be applied
	                                                             ) {
		return detail::lexical_casted_range<T, ForwardRng>( prm_range );
	}

	/// \brief Const range overload of operator| for lexical_casted range adaptor
	template <typename T, typename ForwardRng>
	inline detail::lexical_casted_range<T, const ForwardRng> operator|(const ForwardRng        &prm_range, ///< The range to which the lexical_casted adaptor should be applied
	                                                                   const lexical_casted<T> &           ///< A dummy lexical_casted parameter for determining which adaptor should be applied
	                                                                   ) {
		return detail::lexical_casted_range<T, const ForwardRng>( prm_range );
	}

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_LEXICAL_CASTED_HPP
