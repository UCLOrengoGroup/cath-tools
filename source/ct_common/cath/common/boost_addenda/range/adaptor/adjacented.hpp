/// \file
/// \brief The adjacented class header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_ADJACENTED_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_ADJACENTED_HPP

#include "cath/common/boost_addenda/range/adaptor/range/adjacented_range.hpp"

namespace cath::common {
	namespace detail {

		struct adjacented_forwarder {};

		/// \brief Non-const range overload of operator| for adjacented range adaptor
		template <typename ForwardRng>
		inline adjacented_range<ForwardRng> operator|(ForwardRng &prm_range, ///< The range to which the adjacented adaptor should be applied
		                                              adjacented_forwarder   ///< An equal_grouped_holder parameter for holding the parameters (and for determining which adaptor should be applied)
		                                              ) {
			return adjacented_range<ForwardRng>( prm_range );
		}

		/// \brief Const range overload of operator| for adjacented range adaptor
		template <typename ForwardRng>
		inline adjacented_range<const ForwardRng> operator|(const ForwardRng &prm_range, ///< The range to which the adjacented adaptor should be applied
		                                                    adjacented_forwarder         ///< An equal_grouped_holder parameter for holding the parameters (and for determining which adaptor should be applied)
		                                                    ) {
			return adjacented_range<const ForwardRng>( prm_range );
		}
	} // namespace detail

	inline constexpr detail::adjacented_forwarder adjacented{};

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_ADJACENTED_HPP
