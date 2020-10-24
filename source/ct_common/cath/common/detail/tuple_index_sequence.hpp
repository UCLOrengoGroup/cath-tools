/// \file
/// \brief The tuple_index_sequence header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_DETAIL_TUPLE_INDEX_SEQUENCE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_DETAIL_TUPLE_INDEX_SEQUENCE_HPP

#include <tuple>
#include <type_traits>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Helper to get an index_sequence corresponding to the indices of a tuple type
			template <typename Tpl>
			using tuple_index_sequence = std::make_index_sequence< std::tuple_size< std::decay_t< Tpl > >::value >;

		} // namespace detail
	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_DETAIL_TUPLE_INDEX_SEQUENCE_HPP
