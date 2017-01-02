/// \file
/// \brief The equal_grouped class header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_EQUAL_GROUPED_H
#define _CATH_TOOLS_SOURCE_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_EQUAL_GROUPED_H

#include "common/boost_addenda/range/adaptor/detail/equal_grouped_holder.hpp"
#include "common/boost_addenda/range/adaptor/detail/gen_forwarder.hpp"
#include "common/boost_addenda/range/adaptor/range/equal_grouped_range.hpp"
#include "common/detail/make_static_const.hpp"

namespace cath {
	namespace common {

		/// \brief TODOCUMENT
		class equal_grouped_forwarder final {
		public:
			template <typename FN>
			detail::equal_grouped_holder<FN> operator()(FN) const;

			detail::equal_grouped_holder<std::nullptr_t> operator()() const;
		};

		/// \brief TODOCUMENT
		template <typename FN>
		detail::equal_grouped_holder<FN> equal_grouped_forwarder::operator()(FN arg_unequal_function
		                                                                     ) const {
			return { arg_unequal_function };
		}

		/// \brief TODOCUMENT
		inline detail::equal_grouped_holder<std::nullptr_t> equal_grouped_forwarder::operator()() const {
			return { nullptr };
		}

		namespace detail {

			/// \brief Non-const range overload of operator| for equal_grouped range adaptor
			template <typename RNG, typename FN>
			inline equal_grouped_range<RNG> operator|(RNG                            &arg_range, ///< The range to which the equal_grouped adaptor should be applied
			                                          const equal_grouped_holder<FN> &arg_holder ///< An equal_grouped_holder parameter for holding the parameters (and for determining which adaptor should be applied)
			                                          ) {
				return {
					arg_range,
					arg_holder.template get_function<RNG>()
				};
			}

			/// \brief Const range overload of operator| for equal_grouped range adaptor
			template <typename RNG, typename FN>
			inline equal_grouped_range<const RNG> operator|(const RNG                      &arg_range, ///< The range to which the equal_grouped adaptor should be applied
			                                                const equal_grouped_holder<FN> &arg_holder ///< An equal_grouped_holder parameter for holding the parameters (and for determining which adaptor should be applied)
			                                                ) {
				return {
					arg_range,
					arg_holder.template get_function<const RNG>()
				};
			}
		} // namespace detail
	} // namespace common
} // namespace cath

namespace cath {
	namespace common {

		MAKE_STATIC_CONST( equal_grouped_forwarder, equal_grouped )

	} // namespace common
} // namespace cath

#endif
