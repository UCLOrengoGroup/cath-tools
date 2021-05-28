/// \file
/// \brief The tuple_increment header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_MINS_MAXS_ELEMENT_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_MINS_MAXS_ELEMENT_HPP

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/detail/tuple_index_sequence.hpp"
#include "cath/common/type_traits/is_tuple.hpp"

#include <cassert>
#include <tuple>

namespace cath::common {
	namespace detail {

		/// \brief Update a min and max with a new value
		template <typename T>
		inline int update_min_max_with_value_impl(T       &prm_min,  ///< The min value
		                                          T       &prm_max,  ///< The max value
		                                          const T &prm_value ///< The new value
		                                          ) {
			if ( prm_value < prm_min ) { prm_min = prm_value; }
			if ( prm_value > prm_max ) { prm_max = prm_value; }
			return 0;
		}

		/// \brief Implementation for update_mins_maxs_with_value
		template <typename Tpl, size_t... Index>
		inline void update_mins_maxs_with_value_impl(Tpl       &prm_mins,          ///< The min tuple
		                                             Tpl       &prm_maxs,          ///< The max tuple
		                                             const Tpl &prm_value,         ///< The new tuple
		                                             std::index_sequence<Index...> ///< An index_sequence matching the indices of Tpl
		                                             ) {
			[[maybe_unused]] auto dummy_list = {
				update_min_max_with_value_impl(
					std::get<Index>( prm_mins  ),
					std::get<Index>( prm_maxs  ),
					std::get<Index>( prm_value )
				)...
			};
		}

		/// \brief Element-wise update a min tuple and a max tuple with a new tuple
		template <typename Tpl>
		inline void update_mins_maxs_with_value(Tpl       &prm_mins, ///< The min tuple
		                                        Tpl       &prm_maxs, ///< The max tuple
		                                        const Tpl &prm_value ///< The new tuple
		                                        ) {
			return update_mins_maxs_with_value_impl(
				prm_mins,
				prm_maxs,
				prm_value,
				tuple_index_sequence<Tpl>{}
			);
		}

		/// \brief Function object to find the element-wise mins and maxs tuples of a range of tuples
		struct tuple_mins_maxs_element_fn final {

			/// \brief Find the element-wise mins and maxs tuples of a range of tuples
			template <typename Rng>
			inline auto operator()(Rng &&prm_rng ///< A range of tuples
			                       ) const {
				using value_t = common::range_value_t<Rng>;

				static_assert( is_tuple_v< value_t >, "tuple_mins_maxs_element requires a range of tuples" );

				auto       begin_itr     = ::std::cbegin( prm_rng );
				const auto end_itr       = ::std::cend  ( prm_rng );

				assert( begin_itr != end_itr ); // prm_rng mustn't be empty

				// Set the mins and maxs to the element
				value_t mins = *begin_itr;
				value_t maxs = mins;

				// Update the mins and maxs over the rest of the range
				++begin_itr;
				std::for_each(
					begin_itr,
					end_itr,
					[&] (const value_t &x) { update_mins_maxs_with_value( mins, maxs, x ); }
				);
				return make_pair( mins, maxs );
			}

			tuple_mins_maxs_element_fn()                                   = delete;
			tuple_mins_maxs_element_fn(const tuple_mins_maxs_element_fn &) = delete;
			void operator=(const tuple_mins_maxs_element_fn &)             = delete;
		};

	} // namespace detail

	inline constexpr detail::tuple_mins_maxs_element_fn tuple_mins_maxs_element{};

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_MINS_MAXS_ELEMENT_HPP
