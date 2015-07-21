/// \file
/// \brief The res_pair_keyer class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef RES_PAIR_KEYER_H_INCLUDED
#define RES_PAIR_KEYER_H_INCLUDED

#include "scan/res_pair_keyer/detail/res_pair_keyer_helper.h"
#include "scan/res_pair_keyer/detail/res_pair_keyer_io.h"

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
			using keyer_part_tuple_type = std::tuple<KPs...>;

			/// \brief TODOCUMENT
			keyer_part_tuple_type keyer_parts;

		public:
			/// \brief TODOCUMENT
			using key_tuple_type = detail::key_tuple_t<keyer_part_tuple_type>;

			/// \brief TODOCUMENT
			using key_ranges_tuple_type = detail::key_ranges_tuple_t<keyer_part_tuple_type>;

			res_pair_keyer(const KPs &...);
			// res_pair_keyer(const keyer_part_tuple_type &);

			key_tuple_type make_key(const detail::multi_struc_res_rep_pair &) const;

			key_ranges_tuple_type make_close_keys(const detail::multi_struc_res_rep_pair &,
			                                      const quad_criteria &) const;
			std::string parts_names() const;
		};

		/// \brief TODOCUMENT
		template <typename... KPs>
		inline res_pair_keyer<KPs...>::res_pair_keyer(const KPs &... arg_keyer_parts ///< TODOCUMENT
		                                              ) : keyer_parts( arg_keyer_parts... ) {
		}

		// /// \brief TODOCUMENT
		// template <typename... KPs>
		// inline res_pair_keyer<KPs...>::res_pair_keyer(const keyer_part_tuple_type &arg_keyer_parts ///< TODOCUMENT
		//                                               ) : keyer_parts( arg_keyer_parts ) {
		// }

		/// \brief TODOCUMENT
		template <typename... KPs>
		inline auto res_pair_keyer<KPs...>::make_key(const detail::multi_struc_res_rep_pair &arg_res_pair ///< TODOCUMENT
		                                             ) const -> key_tuple_type {
			return detail::make_key( keyer_parts, arg_res_pair );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		inline auto res_pair_keyer<KPs...>::make_close_keys(const detail::multi_struc_res_rep_pair &arg_res_pair, ///< TODOCUMENT
		                                                    const quad_criteria                    &arg_criteria  ///< The criteria defining what is considered a match
		                                                    ) const -> key_ranges_tuple_type {
			return detail::make_close_keys( keyer_parts, arg_res_pair, arg_criteria );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		std::string res_pair_keyer<KPs...>::parts_names() const {
			return "res_pair_keyer[ " + detail::output_keyer_parts( keyer_parts ) + " ]";
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		res_pair_keyer<KPs...> make_res_pair_keyer(KPs ... arg_keyer_parts ///< TODOCUMENT
		                                           ) {
			return { arg_keyer_parts... };
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		std::ostream & operator<<(std::ostream                 &arg_os,   ///< TODOCUMENT
		                          const res_pair_keyer<KPs...> &arg_keyer ///< TODOCUMENT
		                          ) {
			arg_os << "res_pair_keyer[ " << arg_keyer.parts_names() << " ]";
			return arg_os;
		}

//		template <typename T> class TD;
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

	}
}

#endif
