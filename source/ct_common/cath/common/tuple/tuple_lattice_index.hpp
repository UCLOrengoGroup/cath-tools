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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_LATTICE_INDEX_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_LATTICE_INDEX_HPP

#include <cstddef>
#include <tuple>

#include "cath/common/detail/tuple_index_sequence.hpp"
#include "cath/common/type_traits.hpp"
#include "cath/common/type_traits/is_tuple.hpp"

namespace cath::common {
	namespace detail {

		/// \brief Implementation of tuple_lattice_index
		template <size_t Idx>
		struct tuple_lattice_index_impl final {

			/// \brief Function for implementation of tuple_lattice_index
			template <typename Tpl>
			static constexpr auto fn(const Tpl &prm_indices, ///< The indices of the cell
			                         const Tpl &prm_sizes    ///< The dimensions of the lattice
			                         ) {
				return
					tuple_lattice_index_impl<Idx - 1>::fn( prm_indices, prm_sizes )
					* std::get<Idx>( prm_sizes   )
					+ std::get<Idx>( prm_indices );
			}
		};

		/// \brief Specialisation of implementation of tuple_lattice_index
		template <>
		struct tuple_lattice_index_impl<0> final {

			/// \brief Function for specialisation of implementation of tuple_lattice_index
			template <typename Tpl>
			static constexpr auto fn(const Tpl &prm_indices,  ///< The indices of the cell
			                         const Tpl &/*prm_sizes*/ ///< The dimensions of the lattice
			                         ) {
				return std::get<0>( prm_indices );
			}
		};

		/// \brief Function object to find the overall index of the cell with the specified indices
		///        in a lattice of the specified dimensions
		struct tuple_lattice_index_fn final {

			/// \brief Find the overall index of the cell with the specified indices
			///        in a lattice of the specified dimensions
			template <typename Tpl, typename = std::enable_if< is_tuple_v< Tpl > > >
			constexpr auto operator()(const Tpl &prm_indices, ///< The indices of the cell
			                          const Tpl &prm_sizes    ///< The dimensions of the lattice
			                          ) const {
				static_assert( std::tuple_size_v< common::remove_cvref_t< Tpl > > > 0, "Can't use tuple_lattice_index on tuple with no elements" );
				return tuple_lattice_index_impl<std::tuple_size_v< common::remove_cvref_t< Tpl > > - 1>::fn(
					prm_indices,
					prm_sizes
				);
			}

			tuple_lattice_index_fn()                               = delete;
			tuple_lattice_index_fn(const tuple_lattice_index_fn &) = delete;
			void operator=(const tuple_lattice_index_fn &)         = delete;
		};

	} // namespace detail

	inline constexpr detail::tuple_lattice_index_fn tuple_lattice_index{};

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_LATTICE_INDEX_HPP
