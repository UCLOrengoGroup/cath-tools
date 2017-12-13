/// \file
/// \brief The res_pair_phi_psi_angle_keyer_part class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_DETAIL_RES_PAIR_PHI_PSI_ANGLE_KEYER_PART_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_DETAIL_RES_PAIR_PHI_PSI_ANGLE_KEYER_PART_HPP

#include <boost/lexical_cast.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/join.hpp>

#include "scan/detail/scan_type_aliases.hpp"
#include "structure/geometry/angle.hpp"

#include <string>

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Template for keyer_part that indexes the from/to phi/psi angles
			template <typename F>
			class res_pair_phi_psi_angle_keyer_part final {
			public:
				/// \brief TODOCUMENT
				using value_t           = angle_type;

				/// \brief TODOCUMENT
				using cell_index_t      = key_angle_index_type;

				/// \brief TODOCUMENT
				using cell_index_list_t = boost::range::joined_range<const boost::integer_range<cell_index_t>,
				                                                     const boost::integer_range<cell_index_t>>;

				/// \brief TODOCUMENT
				using search_radius_t   = value_t;

			private:
				/// \brief The cell width to use this dimension
				value_t cell_width;

			public:
				/// \brief Ctor from cell_width
				explicit res_pair_phi_psi_angle_keyer_part(const value_t &arg_cell_width ///< The cell width to use in keying this part
				                                           ) : cell_width( F::sanity_check_cell_width( arg_cell_width ) ) {
				}

				/// \brief Get a short name that describes this key part
				std::string get_name() const {
					return F::get_name() + "(" + boost::lexical_cast<std::string>( cell_width ) + "]";
				}

				/// \brief Extract the relevant value from the specified res_pair
				value_t get_value(const multi_struc_res_rep_pair &arg_res_pair ///< The res_pair from which the relevant value should be extracted
				                     ) const {
					return F::get_value( arg_res_pair );
				}

				/// \brief Extract the search radius from the specified quad_criteria
				search_radius_t get_search_radius(const quad_criteria &arg_criteria   ///< The criteria defining what is considered a match
				                                  ) const {
					return F::get_search_radius( arg_criteria  );
				}

				/// \brief Generate the key part for the specified value
				cell_index_t key_part(const value_t &arg_value ///< The value for which the key_part should be extracted
				                      ) const {
					return static_cast<cell_index_t>( std::floor( arg_value / cell_width ) );
				}

				/// \brief Generate a list of all key parts for all conceivable res_pairs that would match the specified value
				///        within the specified search radius
				cell_index_list_t close_key_parts(const value_t         &arg_value,        ///< The res_pair whose matches' key parts should be generated
				                                  const search_radius_t &arg_search_radius ///< The search radius defining what is considered a match
				                                  ) const {
#ifndef NDEBUG
					// In debug mode, sanity check the inputs
					/// \todo Create a `bool is_shifted(const angle &, ...)` helper function for angle and use it here
					if ( arg_value < geom::zero_angle<angle_base_type>() || arg_value >= geom::one_revolution<angle_base_type>() ) {
						BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot extract key_part for a value that isn't in ( 0, 2pi ]"));
					}
					if ( arg_search_radius <= geom::zero_angle<angle_base_type>() || arg_search_radius >= geom::half_revolution<angle_base_type>() ) {
						BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot search with an angle radius that isn't in ( 0, pi )"));
					}
#endif
					// Calculate:
					//  * the shifted start and stop
					//  * whether this start and stop imply a wrap past 0 degrees
					//  * the begin and one-past-the-end indices
					const auto start           = ( arg_value - arg_search_radius ).quick_shift();
					const auto stop            = ( arg_value + arg_search_radius ).quick_shift();
					const auto wraps           = ( start > stop );
					const cell_index_t begin   = key_part( start );
					const cell_index_t end     = static_cast<cell_index_t>( key_part( stop  ) + 1 );
					const cell_index_t end_all = static_cast<cell_index_t>( std::ceil( geom::one_revolution<angle_base_type>() / cell_width ) );

					// Create a range that will cover the correct cells by joining two integer_ranges.
					// Where the relevant cells wrap past 0 degrees, this looks like [ begin, end_all ) U [ 0, end )
					// Where the relevant cells don't wrap,          this looks like [ begin, end     ) U [ 0, 0   )
					return boost::range::join(
						boost::irange<cell_index_t>( begin, wraps ? end_all : end  ),
						boost::irange<cell_index_t>( 0,     wraps ? end     : 0    )
					);
				}
			};

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
