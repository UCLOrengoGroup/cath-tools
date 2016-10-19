/// \file
/// \brief The hash_tupleheader

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_DETAIL_SCAN_INDEX_STORE_DETAIL_HASH_TUPLE_H
#define _CATH_TOOLS_SOURCE_SCAN_DETAIL_SCAN_INDEX_STORE_DETAIL_HASH_TUPLE_H

#include "scan/detail/res_pair_dirn/res_pair_dirn.h"

#include <functional>

namespace cath {
	namespace scan {
		namespace detail {

			/// The below code attempts to build a hash function for a tuple of standard types
			/// using the code from http://stackoverflow.com/a/7115547 (Leo Goodstadt after Matthieu M.)
			namespace hash_tuple {

				/// \brief Hash class that uses the standard hash for non-tuple types
				template <typename T>
				struct hash {
					size_t operator()(const T &arg_value
					                  ) const {
						return std::hash<T>()( arg_value );
					}
				};

				/// \brief Hash class for res_pair_dirn
				template <>
				struct hash<res_pair_dirn> {
					size_t operator()(const res_pair_dirn &arg_value ///< The res_pair_dirn value
					                  ) const {
						using res_pair_dirn_ut = typename std::underlying_type<res_pair_dirn>::type;
						return std::hash<res_pair_dirn_ut>()( static_cast<res_pair_dirn_ut>( arg_value ) );
					}
				};

				namespace detail {

					/// \brief Method of hashing a new value into existing seed
					///        (using hash_tuple::hash which uses std::hash<>)
					template <class T>
					inline void hash_combine(std::size_t &arg_seed, ///< The seed to update to reflect the addition of arg_value
					                         const T     &arg_value ///< The new value
					                         ) {
						arg_seed ^= hash_tuple::hash<T>()( arg_value ) + 0x9e3779b9 + ( arg_seed << 6 ) + ( arg_seed >> 2 );
					}

					/// \brief Recursive template for hashing tuples
					template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
					struct HashValueImpl final {

						/// \brief Apply the hash method to the existing seed and remaining value
						static void apply(size_t      &seed,
						                  const Tuple &arg_value ///< The new value
						                  ) {
							HashValueImpl<Tuple, Index - 1>::apply( seed, arg_value );
							hash_combine( seed, std::get<Index>( arg_value ) );
						}
					};

					/// \brief End of recursion for hashing tuples
					template <class Tuple>
					struct HashValueImpl<Tuple, 0> {

						/// \brief Apply the has method to the existing seed and the single remaining value of the tuple
						static void apply(size_t      &seed,
						                  const Tuple &arg_value ///< The new value
						                  ) {
							hash_combine( seed, std::get<0>( arg_value ) );
						}
					};
				} // namespace detail

				/// \brief Specialisation of the hash class for tuples
				template <typename ... T>
				struct hash<std::tuple<T...>> {

					/// \brief Hash function that uses the above detail::HashValueImpl code
					size_t operator()(const std::tuple<T...> &arg_value ///< The value to hash
					                  ) const {
						size_t seed = 0;
						detail::HashValueImpl<std::tuple<T...> >::apply( seed, arg_value );
						return seed;
					}
				};

			} // namespace hash_tuple
		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
