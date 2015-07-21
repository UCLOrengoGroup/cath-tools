/// \file
/// \brief The limited class header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Binaries project and then tweaked, eg namespaced in cath)
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

#ifndef LIMITED_H_INCLUDED
#define LIMITED_H_INCLUDED

#include "common/boost_addenda/range/adaptor/detail/gen_forwarder.h"
#include "common/boost_addenda/range/adaptor/detail/limited_holder.h"
#include "common/boost_addenda/range/adaptor/range/limited_range.h"

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Non-const range overload of operator| for limited range adaptor
			template <typename RNG>
			inline limited_range<RNG> operator|(RNG                  &arg_range, ///< The range to which the limited adaptor should be applied
			                                    const limited_holder &arg_holder ///< An limited_holder parameter for holding the parameters (and for determining which adaptor should be applied)
			                                    ) {
				return {
					arg_range,
					arg_holder.get_max_num_elements()
				};
			}

			/// \brief Const range overload of operator| for limited range adaptor
			template <typename RNG>
			inline limited_range<const RNG> operator|(const RNG            &arg_range, ///< The range to which the limited adaptor should be applied
			                                          const limited_holder &arg_holder ///< An limited_holder parameter for holding the parameters (and for determining which adaptor should be applied)
			                                          ) {
				return {
					arg_range,
					arg_holder.get_max_num_elements()
				};
			}
		}
	}
}

namespace cath {
	namespace common {
		namespace {

			/// \brief Following Boost Range's adaptor implementations, create a static forwarder
			///        whose function operator constructs an limited_holder
			static cath::common::detail::gen_forwarder<detail::limited_holder> limited
				= cath::common::detail::gen_forwarder<detail::limited_holder>();
		}
	}
}

#endif
