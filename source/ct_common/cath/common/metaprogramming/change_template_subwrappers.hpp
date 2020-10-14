/// \file
/// \brief The change_template_subwrappers header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_METAPROGRAMMING_CHANGE_TEMPLATE_SUBWRAPPERS_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_METAPROGRAMMING_CHANGE_TEMPLATE_SUBWRAPPERS_HPP

#include "cath/common/metaprogramming/change_template_wrapper.hpp"

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Primary template for implementation of change_template_subwrappers_t
			template <typename T,
			          template <typename...> class NewWrppr>
			struct change_template_subwrappers final {};

			/// \brief Specialisation of template for implementation of change_template_subwrappers_t
			template <template <typename...> class OuterWrppr,
			          typename... Prms,
			          template <typename...> class NewWrppr>
			struct change_template_subwrappers<OuterWrppr<Prms...>, NewWrppr> final {
				using type = OuterWrppr< change_template_wrapper_t<Prms, NewWrppr>... >;
			};

		} // namespace detail

	/// \brief The type created by wrapping the specified type's template parameters in the template wrapper
	///
	/// Eg Given `tuple< tuple<int>, list<double> >` and `vector` this returns `tuple< vector<int>, vector<double> >`
	template <typename T, template <typename...> class Tmplt>
	using change_template_subwrappers_t = typename detail::change_template_subwrappers<T, Tmplt>::type;

	} // namespace common
} // namespace cath

#endif
