/// \file
/// \brief The is_tuple header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_TUPLE_TUPLE_LATTICE_INDEX_H
#define _CATH_TOOLS_SOURCE_COMMON_TUPLE_TUPLE_LATTICE_INDEX_H

#include "common/detail/make_static_const.hpp"
#include "common/detail/tuple_index_sequence.hpp"
#include "common/tuple/is_tuple.hpp"

#include <tuple>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Implementation of tuple_lattice_index
			template <size_t Idx>
			struct tuple_lattice_index_impl final {

				/// \brief Function for implementation of tuple_lattice_index
				template <typename Tpl>
				static constexpr auto fn(const Tpl &arg_indices, ///< The indices of the cell
				                         const Tpl &arg_sizes    ///< The dimensions of the lattice
				                         ) {
					return
						tuple_lattice_index_impl<Idx - 1>::fn( arg_indices, arg_sizes )
						* std::get<Idx>( arg_sizes   )
						+ std::get<Idx>( arg_indices );
				}
			};

			/// \brief Specialisation of implementation of tuple_lattice_index
			template <>
			struct tuple_lattice_index_impl<0> final {

				/// \brief Function for specialisation of implementation of tuple_lattice_index
				template <typename Tpl>
				static constexpr auto fn(const Tpl &arg_indices,  ///< The indices of the cell
				                         const Tpl &/*arg_sizes*/ ///< The dimensions of the lattice
				                         ) {
					return std::get<0>( arg_indices );
				}
			};

			/// \brief Function object to find the overall index of the cell with the specified indices
			///        in a lattice of the specified dimensions
			struct tuple_lattice_index_fn final {

				/// \brief Find the overall index of the cell with the specified indices
				///        in a lattice of the specified dimensions
				template <typename Tpl, typename = std::enable_if< is_tuple< Tpl >::value > >
				constexpr auto operator()(const Tpl &arg_indices, ///< The indices of the cell
				                          const Tpl &arg_sizes    ///< The dimensions of the lattice
				                          ) const {
					constexpr size_t the_tuple_size = std::tuple_size< std::decay_t< Tpl > >::value;
					static_assert( the_tuple_size > 0, "Can't use tuple_lattice_index on tuple with no elements" );
					return tuple_lattice_index_impl<std::tuple_size< std::decay_t< Tpl > >::value - 1>::fn(
						arg_indices,
						arg_sizes
					);
				}

				tuple_lattice_index_fn()                               = delete;
				tuple_lattice_index_fn(const tuple_lattice_index_fn &) = delete;
				void operator=(const tuple_lattice_index_fn &)         = delete;
			};

		} // namespace detail

		MAKE_STATIC_CONST( detail::tuple_lattice_index_fn, tuple_lattice_index )

	} // namespace common
} // namespace cath

#endif
