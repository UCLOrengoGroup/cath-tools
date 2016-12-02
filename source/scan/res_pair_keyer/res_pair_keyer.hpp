/// \file
/// \brief The res_pair_keyer class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_H
#define _CATH_TOOLS_SOURCE_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_H

#include "scan/res_pair_keyer/detail/res_pair_keyer_helper.hpp"
#include "scan/res_pair_keyer/detail/res_pair_keyer_io.hpp"

#include <ostream>
#include <tuple>

namespace cath { namespace scan { class quad_criteria; } }
namespace cath { namespace scan { namespace detail { class multi_struc_res_rep_pair; } } }

namespace cath {
	namespace scan {

		/// \brief Extract the multi_struc_res_rep_pair
		///
		/// \tparam TPL must be a std::tuple of res_pair_keyer_parts
		///
		/// \todo Create and check against a res_pair_keyer_part concept
		template <typename... KPs>
		class res_pair_keyer final {
		private:
			// /// \brief Compile-time check that TPL's a std::tuple<> of zero or more types
			// static_assert( common::is_tuple<TPL>::value, "ERROR: res_pair_keyer must be instantiated with a type that's a std::tuple" );

			/// \brief TODOCUMENT
			using keyer_part_tuple = std::tuple<KPs...>;

			/// \brief TODOCUMENT
			keyer_part_tuple keyer_parts;

		public:
			/// \brief TODOCUMENT
			using key_value_tuple_type = detail::key_value_tuple_t<KPs...>;

			/// \brief TODOCUMENT
			using key_index_tuple_type = detail::key_index_tuple_t<KPs...>;

			/// \brief TODOCUMENT
			using key_ranges_tuple_type = detail::key_ranges_tuple_t<KPs...>;

			explicit constexpr res_pair_keyer(const KPs &...);
			// res_pair_keyer(const keyer_part_tuple &);

			template <typename Data>
			constexpr key_value_tuple_type make_value(Data &&) const;

			template <typename Data>
			constexpr key_index_tuple_type make_key(Data &&) const;

			template <typename Data, typename Crit>
			constexpr key_ranges_tuple_type make_close_keys(Data &&,
			                                                Crit &&) const;

			std::string parts_names() const;
		};

		/// \brief TODOCUMENT
		template <typename... KPs>
		inline constexpr res_pair_keyer<KPs...>::res_pair_keyer(const KPs &... arg_keyer_parts ///< TODOCUMENT
		                                                        ) : keyer_parts( arg_keyer_parts... ) {
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		template <typename Data>
		inline constexpr auto res_pair_keyer<KPs...>::make_value(Data &&arg_data ///< TODOCUMENT
		                                                         ) const -> key_value_tuple_type {
			return detail::make_value( keyer_parts, std::forward<Data>( arg_data ) );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		template <typename Data>
		inline constexpr auto res_pair_keyer<KPs...>::make_key(Data &&arg_data ///< TODOCUMENT
		                                                       ) const -> key_index_tuple_type {
			return detail::make_key( keyer_parts, std::forward<Data>( arg_data ) );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		template <typename Data, typename Crit>
		inline constexpr auto res_pair_keyer<KPs...>::make_close_keys(Data &&arg_data, ///< TODOCUMENT
		                                                              Crit &&arg_criteria  ///< The criteria defining what is considered a match
		                                                              ) const -> key_ranges_tuple_type {
			return detail::make_close_keys(
				keyer_parts,
				std::forward<Data>( arg_data ),
				std::forward<Crit>( arg_criteria )
			);
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		std::string res_pair_keyer<KPs...>::parts_names() const {
			return "res_pair_keyer[ " + detail::output_keyer_parts( keyer_parts ) + " ]";
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		constexpr res_pair_keyer<KPs...> make_res_pair_keyer(KPs ... arg_keyer_parts ///< TODOCUMENT
		                                                     ) {
			/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
			return res_pair_keyer<KPs...>{ arg_keyer_parts... };
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		std::ostream & operator<<(std::ostream                 &arg_os,   ///< TODOCUMENT
		                          const res_pair_keyer<KPs...> &arg_keyer ///< TODOCUMENT
		                          ) {
			arg_os << "res_pair_keyer[ " << arg_keyer.parts_names() << " ]";
			return arg_os;
		}

//		auto make_example() {
//			auto the_example = make_res_pair_keyer(
//				res_pair_from_phi_keyer_part  ( geom::make_angle_from_degrees<detail::angle_base_type>( 67.5 ) ),
//				res_pair_view_x_keyer_part    ( 2.0 ),
////				res_pair_orient_keyer_part( ),
//				res_pair_index_dirn_keyer_part( )
//
//			);
////			TD< decltype( the_example ) > bob;
////			TD< decltype( the_example.make_key( std::declval<multi_struc_res_rep_pair>(), std::declval<quad_criteria>() ) ) > fred;
//			TD< decltype(the_example)::key_ranges_tuple_type > fred;
//			return the_example;
//		}

	} // namespace scan
} // namespace cath

#endif
