/// \file
/// \brief The limited class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_LIMITED_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_LIMITED_HPP

#include <boost/core/ignore_unused.hpp>

#include "cath/common/boost_addenda/range/adaptor/detail/gen_forwarder.hpp"
#include "cath/common/boost_addenda/range/adaptor/detail/limited_holder.hpp"
#include "cath/common/boost_addenda/range/adaptor/range/limited_range.hpp"

namespace cath::common {
	namespace detail {

		/// \brief Non-const range overload of operator| for limited range adaptor
		template <typename RNG>
		inline limited_range<RNG> operator|(RNG                  &prm_range, ///< The range to which the limited adaptor should be applied
		                                    const limited_holder &prm_holder ///< An limited_holder parameter for holding the parameters (and for determining which adaptor should be applied)
		                                    ) {
			return {
				prm_range,
				prm_holder.get_max_num_elements()
			};
		}

		/// \brief Const range overload of operator| for limited range adaptor
		template <typename RNG>
		inline limited_range<const RNG> operator|(const RNG            &prm_range, ///< The range to which the limited adaptor should be applied
		                                          const limited_holder &prm_holder ///< An limited_holder parameter for holding the parameters (and for determining which adaptor should be applied)
		                                          ) {
			return {
				prm_range,
				prm_holder.get_max_num_elements()
			};
		}
	} // namespace detail

	inline constexpr detail::gen_forwarder<detail::limited_holder> limited{};

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_LIMITED_HPP
