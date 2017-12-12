/// \file
/// \brief The combine_params_lists_with_template_list header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_METAPROGRAMMING_COMBINE_PARAMS_LISTS_WITH_TEMPLATE_LIST_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_METAPROGRAMMING_COMBINE_PARAMS_LISTS_WITH_TEMPLATE_LIST_H

#include "common/metaprogramming/append_template_params_into_first_wrapper.hpp"
#include "common/metaprogramming/change_template_subwrappers.hpp"

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Primary template for implementation of combine_params_lists_with_template_list_t
			template <typename T,
			          typename U>
			struct combine_params_lists_with_template_list final {};

			/// \brief Specialisation of template for implementation of combine_params_lists_with_template_list_t
			template <typename U,
			          template <template <typename...> class...> class TmpltWrp,
			          template <typename...> class... Tmplts
			          >
			struct combine_params_lists_with_template_list<U,
			                                               TmpltWrp<Tmplts...>
			                                               > final {
				using type = append_template_params_into_first_wrapper_t<
					change_template_subwrappers_t<U, Tmplts>...
				>;
			};

		} // namespace detail

	/// \brief The type created by rewrapping each of the typelists types in the specified typelist with
	///        each of the templates in the specified template list and then rewrapping the (flattened) list
	///        in the typelist wrapper used for the specified typelist.
	///
	/// Eg Given `tuple< tuple<int>, tuple<double> >` and `template_list<vector, list>` this returns 
	/// `tuple< vector<int>, vector<double>, list<int>, list<double> >`
	template <typename T,
	          typename U>
	using combine_params_lists_with_template_list_t = typename detail::combine_params_lists_with_template_list<T, U>::type;

	} // namespace common
} // namespace cath

#endif
