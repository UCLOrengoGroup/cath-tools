/// \file
/// \brief The res_pair_keyer helper header

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
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

#ifndef _CATH_TOOLS_SOURCE_SCAN_RES_PAIR_KEYER_DETAIL_RES_PAIR_KEYER_HELPER_H
#define _CATH_TOOLS_SOURCE_SCAN_RES_PAIR_KEYER_DETAIL_RES_PAIR_KEYER_HELPER_H

#include "common/cpp17/apply.h"

namespace cath { namespace scan { namespace detail { class multi_struc_res_rep_pair; } } }
namespace cath { namespace scan { class quad_criteria; } }

#include <utility>

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Metafunction to extract value_t from a keyer_part
			template <typename T>
			using keyer_part_value_type_t           = typename T::value_t;

			/// \brief Metafunction to extract cell_index_t from a keyer_part
			template <typename T>
			using keyer_part_cell_index_type_t      = typename T::cell_index_t;

			/// \brief Metafunction to extract cell_index_list_t from a keyer_part
			template <typename T>
			using keyer_part_cell_index_list_type_t = typename T::cell_index_list_t;



			/// \brief Use the specified keyer_part to extract a key_part from the specified res_pair
			///
			/// This is a helper function for calling keyer_part's key_part() with a
			/// multi_struc_res_rep_pair (rather than with a value)
			template <typename KP>
			auto key_part(const KP                       &arg_keyer_part, ///< The keyer_part with which to extract the key part from the res_pair
			              const multi_struc_res_rep_pair &arg_res_pair    ///< The res_pair whose key part should be extracted
			              ) {
				return arg_keyer_part.key_part( arg_keyer_part.get_value( arg_res_pair ) );
			}


			/// \brief Use the specified keyer_part to extract a list of all key parts for
			///        all conceivable res_pairs that would match the specified res_pair
			///
			/// This is a helper function for calling keyer_part's close_key_parts() with a
			/// multi_struc_res_rep_pair and quad_criteria (rather than with a value and search_radius)
			template <typename KP>
			auto close_key_parts(const KP                       &arg_keyer_part, ///< The keyer_part with which to extract the list of key parts from the res_pair
			                     const multi_struc_res_rep_pair &arg_res_pair,   ///< The res_pair whose matches' key parts should be generated
			                     const quad_criteria            &arg_criteria    ///< The criteria defining what is considered a match
			                     ) {
				return arg_keyer_part.close_key_parts(
					arg_keyer_part.get_value        ( arg_res_pair ),
					arg_keyer_part.get_search_radius( arg_criteria )
				);
			}


			/// \brief Metafunction for calculating the key tuple type from a tuple of key_part types
			///
			/// This is just a definition of the template to be specialised below
			template <typename T> struct key_tuple;

			/// \brief Metafunction for calculating the key tuple type from a tuple of key_part types
			///
			/// This specialisation of the above template does the real work
			template <typename... Ts>
			struct key_tuple<std::tuple<Ts...>> final {
				using type = std::tuple<keyer_part_cell_index_type_t<Ts>...>;
			};

			/// \brief Convenience type alias for metafunction to calculate the key tuple type from a tuple of key_part types
			template <typename T>
			using key_tuple_t = typename key_tuple<T>::type;



			/// \brief Metafunction for calculating the key tuple type from a tuple of key_part types
			///
			/// This is just a definition of the template to be specialised below
			template <typename T> struct key_ranges_tuple;

			/// \brief Metafunction for calculating the key tuple type from a tuple of key_part types
			///
			/// This specialisation of the above template does the real work
			template <typename... Ts>
			struct key_ranges_tuple<std::tuple<Ts...>> final {
				using type = std::tuple<keyer_part_cell_index_list_type_t<Ts>...>;
			};

			/// \brief Convenience type alias for metafunction to calculate the key tuple type from a tuple of key_part types
			template <typename T>
			using key_ranges_tuple_t = typename key_ranges_tuple<T>::type;



			/// \brief TODOCUMENT
			class keyer_part_maker final {
			private:
				/// \brief TODOCUMENT
				const multi_struc_res_rep_pair &the_res_pair;

			public:
				/// \brief TODOCUMENT
				keyer_part_maker(const multi_struc_res_rep_pair &arg_res_pair ///< TODOCUMENT
				                 ) : the_res_pair ( arg_res_pair ) {
				}

				/// \brief TODOCUMENT
				template <typename... Ts>
				constexpr auto operator()(const Ts &... arg_keyer_parts ///< TODOCUMENT
				                          ) {
					return std::make_tuple( key_part( arg_keyer_parts, the_res_pair )... );
				}
			};

			/// \brief TODOCUMENT
			class keyer_part_range_maker final {
			private:
				/// \brief TODOCUMENT
				const multi_struc_res_rep_pair &the_res_pair;

				/// \brief TODOCUMENT
				const quad_criteria &the_criteria;

			public:
				/// \brief TODOCUMENT
				keyer_part_range_maker(const multi_struc_res_rep_pair &arg_res_pair, ///< TODOCUMENT
				                       const quad_criteria            &arg_criteria  ///< TODOCUMENT
				                       ) : the_res_pair ( arg_res_pair ),
				                           the_criteria ( arg_criteria ) {
				}

				/// \brief TODOCUMENT
				template <typename... Ts>
				constexpr auto operator()(const Ts &... arg_keyer_parts
				                          ) {
					return std::make_tuple( close_key_parts( arg_keyer_parts, the_res_pair, the_criteria )... );
				}
			};

			/// \brief Template function to make a key from TODOCUMENT
			template <typename... Ts>
			auto make_key(const std::tuple<Ts...>        &arg_tuple,   ///< The tuple of keyer_parts to be used to make the key
			              const multi_struc_res_rep_pair &arg_res_pair ///< The res_pair whose matches' key parts should be generated
			              ) {
				return common::apply( keyer_part_maker( arg_res_pair ), arg_tuple );
			}

			/// \brief TODOCUMENT
			template <typename... Ts>
			auto make_close_keys(const std::tuple<Ts...>        &arg_tuple,    ///< TODOCUMENT
			                     const multi_struc_res_rep_pair &arg_res_pair, ///< The res_pair whose matches' key parts should be generated
			                     const quad_criteria            &arg_criteria  ///< The criteria defining what is considered a match
			                     ) {
				return common::apply( keyer_part_range_maker( arg_res_pair, arg_criteria ), arg_tuple );
			}

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
