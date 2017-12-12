/// \file
/// \brief The gen_forwarder class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_DETAIL_GEN_FORWARDER_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_DETAIL_GEN_FORWARDER_H

namespace cath {
	namespace common {
		namespace detail {

			/// \brief gen_forwarder used to implement range adaptors
			///
			/// This is a variadic template version of the forwarder classes used in the Boost Range adaptor implementations
			template <typename Holder>
			struct gen_forwarder final {

				/// \brief Return the result of forwarding the parameters to the relevant Holder ctor
				template <typename... T>
				Holder operator()(T ...t) const {
					return Holder( t... );
				}
			};

		} // namespace detail
	} // namespace common
} // namespace cath

#endif
