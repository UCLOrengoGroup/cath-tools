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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_DETAIL_RES_PAIR_KEYER_HELPER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_DETAIL_RES_PAIR_KEYER_HELPER_HPP

#include "cath/common/cpp17/apply.hpp"
#include "cath/common/tuple/make_tuple_with_skips.hpp"

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


			/// \brief Convenience type alias for metafunction to calculate the key tuple type from key_part types
			template <typename... KPs>
			using key_value_tuple_t = common::tuple_with_skips_t<keyer_part_value_type_t<KPs>...>;

			/// \brief Convenience type alias for metafunction to calculate the key tuple type from key_part types
			template <typename... KPs>
			using key_index_tuple_t = common::tuple_with_skips_t<keyer_part_cell_index_type_t<KPs>...>;

			/// \brief Convenience type alias for metafunction to calculate the key tuple type from key_part types
			template <typename... KPs>
			using key_ranges_tuple_t = common::tuple_with_skips_t<keyer_part_cell_index_list_type_t<KPs>...>;



			/// \brief Helper for make_value(), below, (via make_keyer_parts_value_maker())
			///
			/// \todo Come C++17, this can be replaced by a constexpr lambda - see make_value() notes
			template <typename Data>
			struct keyer_parts_value_maker final {

				/// \brief Const reference to the data to be passed to the keyer_parts
				///
				/// This mustn't attempt to store by rvalue ref and forward because it's
				/// used in multiple calls for the keyer_parts in the tuple
				const Data &data;
				
				/// \brief Ctor from the data to be passed to the keyer_parts
				explicit constexpr keyer_parts_value_maker(const Data &prm_data ///< The data to be passed to the keyer_parts
				                                           ) : data{ prm_data } {
				}
				
				/// \brief Function operator to apply a bunch of keyer_parts to data and return the values in a tuple
				template <typename... KPs>
				constexpr auto operator()(const KPs &...prm_keyer_parts ///< The keyer_parts to be used to get the value parts from the data
				                          ) {
					return common::make_tuple_with_skips( prm_keyer_parts.get_value( data )... );
				}
			};

			/// \brief Factory function for keyer_parts_value_maker to provide type deduction for construction
			///
			/// \todo Come C++17 there are two reasons this is redundant and can be removed:
			///        * class template deduction means the deduction can just be performed by a call to the ctor
			///        * keyer_parts_value_maker can be completely replaced with a constexpr lambda
			template <typename Data>
			constexpr keyer_parts_value_maker<Data> make_keyer_parts_value_maker(const Data &prm_data ///< The data to be passed to the keyer_parts
			                                                                     ) {
				/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
				return keyer_parts_value_maker<Data>{ prm_data };
			}





			/// \brief Helper for make_key(), below, (via make_keyer_parts_key_maker())
			///
			/// \todo Come C++17, this can be replaced by a constexpr lambda - see make_key() notes
			template <typename Data>
			struct keyer_parts_key_maker final {

				/// \brief Const reference to the data to be passed to the keyer_parts
				///
				/// This mustn't attempt to store by rvalue ref and forward because it's
				/// used in multiple calls for the keyer_parts in the tuple
				const Data &data;
				
				/// \brief Ctor from the data to be passed to the keyer_parts
				explicit constexpr keyer_parts_key_maker(const Data &prm_data ///< The data to be passed to the keyer_parts
				                                         ) : data{ prm_data } {
				}
				
				/// \brief Function operator to apply a bunch of keyer_parts to data and return the keys in a tuple
				template <typename... KPs>
				constexpr auto operator()(const KPs &...prm_keyer_parts ///< The keyer_parts to be used to get the key parts from the data
				                          ) {
					return common::make_tuple_with_skips(
						prm_keyer_parts.key_part(
							prm_keyer_parts.get_value( data )
						)...
					);
				}
			};

			/// \brief Factory function for keyer_parts_key_maker to provide type deduction for construction
			///
			/// \todo Come C++17 there are two reasons this is redundant and can be removed:
			///        * class template deduction means the deduction can just be performed by a call to the ctor
			///        * keyer_parts_key_maker can be completely replaced with a constexpr lambda
			template <typename Data>
			constexpr keyer_parts_key_maker<Data> make_keyer_parts_key_maker(const Data &prm_data ///< The data to be passed to the keyer_parts
			                                                                 ) {
				/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
				return keyer_parts_key_maker<Data>{ prm_data };
			}






			/// \brief Helper for make_min_close_key(), below, (via make_keyer_parts_min_close_key_maker())
			///
			/// \todo Come C++17, this can be replaced by a constexpr lambda - see make_min_close_key() notes
			template <typename Data, typename Crit>
			struct keyer_parts_min_close_key_maker final {

				/// \brief Const reference to the data to be passed to the keyer_parts
				///
				/// This mustn't attempt to store by rvalue ref and forward because it's
				/// used in multiple calls for the keyer_parts in the tuple
				const Data &data;

				/// \brief Const reference to the data to be passed to the keyer_parts
				///
				/// This mustn't attempt to store by rvalue ref and forward because it's
				/// used in multiple calls for the keyer_parts in the tuple
				const Crit &crit;
				
				/// \brief Ctor from the data to be passed to the keyer_parts
				constexpr keyer_parts_min_close_key_maker(const Data &prm_data, ///< The data to be passed to the keyer_parts
				                                          const Crit &prm_crit  ///< TODOCUMENT
				                                          ) : data{ prm_data },
				                                              crit{ prm_crit } {
				}
				
				/// \brief Function operator to apply a bunch of keyer_parts to data and return the keys in a tuple
				template <typename... KPs>
				constexpr auto operator()(const KPs &...prm_keyer_parts ///< The keyer_parts to be used to get the key parts from the data
				                          ) {
					return common::make_tuple_with_skips(
						prm_keyer_parts.min_close_key_part(
							prm_keyer_parts.get_value        ( data ),
							prm_keyer_parts.get_search_radius( crit )
						)...
					);
				}
			};

			/// \brief Factory function for keyer_parts_min_close_key_maker to provide type deduction for construction
			///
			/// \todo Come C++17 there are two reasons this is redundant and can be removed:
			///        * class template deduction means the deduction can just be performed by a call to the ctor
			///        * keyer_parts_min_close_key_maker can be completely replaced with a constexpr lambda
			template <typename Data, typename Crit>
			constexpr keyer_parts_min_close_key_maker<Data, Crit> make_keyer_parts_min_close_key_maker(const Data &prm_data, ///< The data to be passed to the keyer_parts
			                                                                                           const Crit &prm_crit  ///< TODOCUMENT
			                                                                                           ) {
				return { prm_data, prm_crit };
			}




			/// \brief Helper for make_max_close_key(), below, (via make_keyer_parts_max_close_key_maker())
			///
			/// \todo Come C++17, this can be replaced by a constexpr lambda - see make_max_close_key() notes
			template <typename Data, typename Crit>
			struct keyer_parts_max_close_key_maker final {

				/// \brief Const reference to the data to be passed to the keyer_parts
				///
				/// This mustn't attempt to store by rvalue ref and forward because it's
				/// used in multiple calls for the keyer_parts in the tuple
				const Data &data;

				/// \brief Const reference to the data to be passed to the keyer_parts
				///
				/// This mustn't attempt to store by rvalue ref and forward because it's
				/// used in multiple calls for the keyer_parts in the tuple
				const Crit &crit;
				
				/// \brief Ctor from the data to be passed to the keyer_parts
				constexpr keyer_parts_max_close_key_maker(const Data &prm_data, ///< The data to be passed to the keyer_parts
				                                          const Crit &prm_crit  ///< TODOCUMENT
				                                          ) : data{ prm_data },
				                                              crit{ prm_crit } {
				}
				
				/// \brief Function operator to apply a bunch of keyer_parts to data and return the keys in a tuple
				template <typename... KPs>
				constexpr auto operator()(const KPs &...prm_keyer_parts ///< The keyer_parts to be used to get the key parts from the data
				                          ) {
					return common::make_tuple_with_skips(
						prm_keyer_parts.max_close_key_part(
							prm_keyer_parts.get_value        ( data ),
							prm_keyer_parts.get_search_radius( crit )
						)...
					);
				}
			};

			/// \brief Factory function for keyer_parts_max_close_key_maker to provide type deduction for construction
			///
			/// \todo Come C++17 there are two reasons this is redundant and can be removed:
			///        * class template deduction means the deduction can just be performed by a call to the ctor
			///        * keyer_parts_max_close_key_maker can be completely replaced with a constexpr lambda
			template <typename Data, typename Crit>
			constexpr keyer_parts_max_close_key_maker<Data, Crit> make_keyer_parts_max_close_key_maker(const Data &prm_data, ///< The data to be passed to the keyer_parts
			                                                                                           const Crit &prm_crit  ///< TODOCUMENT
			                                                                                           ) {
				return { prm_data, prm_crit };
			}







			/// \brief Template function to make a value tuple from a tuple of keyer_parts and an instance of their common data type
			///
			/// \todo Come C++17, remove keyer_parts_value_maker and replace the make_keyer_parts_value_maker() call
			///       below with a constexpr lambda like:
			///     [&] (const auto &...keyer_parts) { return std::make_tuple( keyer_parts.get_value( prm_data )... ); }
			///
			template <typename... KPs, typename Data>
			constexpr decltype(auto) make_value(const std::tuple<KPs...>  &prm_tuple, ///< The tuple of keyer_parts to be used to make the value
			                                    const Data                &prm_data   ///< The data to be passed to the keyer_parts
			                                    ) {
				return common::apply( make_keyer_parts_value_maker( prm_data ), prm_tuple );
			}

			/// \brief Template function to emplace_back in a store the value components from a tuple of keyer_parts and an instance of their common data type
			///
			/// Unlike make_value() and make_key(), this just uses a lambda because the return type's won't permit constexpr anyway.
			///
			/// This mustn't attempt to store prm_data by rvalue ref and forward because it's being
			/// forwarded to multiple calls for the keyer_parts in the tuple
			template <typename Store, typename Key, typename... KPs, typename Data>
			inline void store_emplace_value(Store                     &prm_store, ///< The store in which to emplace_back the value components
			                                const Key                 &prm_key,   ///< The key under which the value should be recorded
			                                const std::tuple<KPs...>  &prm_tuple, ///< The tuple of keyer_parts to apply to be used to make the close_keys
			                                const Data                &prm_data   ///< The data to be passed to the keyer_parts
			                                ) {
				common::apply(
					[&] (const auto &...keyer_parts) {
						prm_store.emplace_back_entry_to_cell(
							prm_key,
							keyer_parts.get_value( prm_data )...
						);
					},
					prm_tuple
				);
			}

			/// \brief Template function to make a key tuple from a tuple of keyer_parts and an instance of their common data type
			///
			/// \todo Come C++17, remove keyer_parts_key_maker and replace the make_keyer_parts_key_maker() call
			///       below with a constexpr lambda like:
			///
			///     [&] (const auto &...keyer_parts) { return std::make_tuple( keyer_parts.key_part( keyer_parts.get_value( prm_data ) )... ); }
			template <typename... KPs, typename Data>
			constexpr decltype(auto) make_key(const std::tuple<KPs...>  &prm_tuple, ///< The tuple of keyer_parts to be used to make the key
			                                  const Data                &prm_data   ///< The data to be passed to the keyer_parts
			                                  ) {
				return common::apply( make_keyer_parts_key_maker( prm_data ), prm_tuple );
			}

			/// \brief Template function to make a close_keys tuple of ranges from a tuple of keyer_parts and an instance of their common data type
			///
			/// Unlike make_value() and make_key(), this just uses a lambda because the return type's won't permit constexpr anyway.
			///
			/// This mustn't attempt to store prm_data/prm_crit by rvalue ref and forward because they're being
			/// forwarded to multiple calls for the keyer_parts in the tuple
			template <typename... KPs, typename Data, typename Crit>
			decltype(auto) make_close_keys(const std::tuple<KPs...>  &prm_tuple, ///< The tuple of keyer_parts to apply to be used to make the close_keys
			                               const Data                &prm_data,  ///< The data to be passed to the keyer_parts
			                               const Crit                &prm_crit   ///< The criteria defining what is considered a match
			                               ) {
				return common::apply(
					[&] (const auto &...keyer_parts) {
						return common::make_tuple_with_skips( keyer_parts.close_key_parts(
							keyer_parts.get_value        ( prm_data ),
							keyer_parts.get_search_radius( prm_crit )
						)... );
					},
					prm_tuple
				);
			}

			/// \brief Template function to make a min-close-key tuple from a tuple of keyer_parts and an instance of their common data type
			///
			/// \todo Come C++17, remove keyer_parts_key_maker and replace the make_keyer_parts_key_maker() call
			///       below with a constexpr lambda like:
			///
			///     [&] (const auto &...keyer_parts) {
			///       return std::make_tuple(
			///         keyer_parts.min_close_key_part(
			///           keyer_parts.get_value        ( prm_data ),
			///           keyer_parts.get_search_radius( prm_crit )
			///         )...
			///       );
			///     }
			template <typename... KPs, typename Data, typename Crit>
			constexpr decltype(auto) make_min_close_key(const std::tuple<KPs...>  &prm_tuple, ///< The tuple of keyer_parts to be used to make the key
			                                            const Data                &prm_data,  ///< The data to be passed to the keyer_parts
			                                            const Crit                &prm_crit   ///< The criteria defining what is considered a match
			                                            ) {
				return common::apply( make_keyer_parts_min_close_key_maker( prm_data, prm_crit ), prm_tuple );
			}

			/// \brief Template function to make a max-close-key tuple from a tuple of keyer_parts and an instance of their common data type
			///
			/// \todo Come C++17, remove keyer_parts_key_maker and replace the make_keyer_parts_key_maker() call
			///       below with a constexpr lambda like:
			///
			///     [&] (const auto &...keyer_parts) {
			///       return std::make_tuple(
			///         keyer_parts.max_close_key_part(
			///           keyer_parts.get_value        ( prm_data ),
			///           keyer_parts.get_search_radius( prm_crit )
			///         )...
			///       );
			///     }
			template <typename... KPs, typename Data, typename Crit>
			constexpr decltype(auto) make_max_close_key(const std::tuple<KPs...>  &prm_tuple, ///< The tuple of keyer_parts to be used to make the key
			                                            const Data                &prm_data,  ///< The data to be passed to the keyer_parts
			                                            const Crit                &prm_crit   ///< The criteria defining what is considered a match
			                                            ) {
				return common::apply( make_keyer_parts_max_close_key_maker( prm_data, prm_crit ), prm_tuple );
			}

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_DETAIL_RES_PAIR_KEYER_HELPER_HPP
