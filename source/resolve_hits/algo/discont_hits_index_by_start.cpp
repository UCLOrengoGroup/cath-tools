/// \file
/// \brief The discont_hits_index_by_start class definitions

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

#include "discont_hits_index_by_start.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/algorithm/upper_bound.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/copy_build.h"
#include "common/algorithm/sort_uniq_copy.h"

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::integer_range;
using boost::irange;
using boost::numeric_cast;
using boost::range::lower_bound;
using boost::range::upper_bound;
using std::make_pair;
using std::string;

/// \brief Calculate the starts and indices of the discontiguous domains in the specified calc_hit_list,
///        ordered by their starts
res_arr_idx_pair_vec discont_hits_index_by_start::calc_disconts(const calc_hit_list &arg_calc_hit_list ///< The calc_hit_list to index
                                                                ) {
	return sort_copy(
		copy_build<res_arr_idx_pair_vec>(
			irange( 0_z, arg_calc_hit_list.size() )
				| filtered   ( [&] (const size_t &x) {
					return arg_calc_hit_list[ x ].is_discontig();
				} )
				| transformed( [&] (const size_t &x) {
					return make_pair( arg_calc_hit_list[ x ].get_start_arrow(), numeric_cast<hitidx_t>( x ) );
				} )
		)
	);
}

/// \brief Ctor from a calc_hit_list
discont_hits_index_by_start::discont_hits_index_by_start(const calc_hit_list &arg_calc_hit_list ///< The calc_hit_list to index
                                                         ) : the_hits ( arg_calc_hit_list                   ),
                                                             disconts ( calc_disconts ( arg_calc_hit_list ) ) {
}

/// \brief Get this index's internal indices for any discontiguous domains that start
///        within the specified region (inclusive)
///
/// The only way these indices should be used is as arguments to get_discont_hit_of_index_index()
integer_range<size_t> discont_hits_index_by_start::get_index_indices_of_disconts_in_range(const res_arrow &arg_start, ///< The start of the region of interest
                                                                                          const res_arrow &arg_stop   ///< The end of the region of interest
                                                                                          ) const {
	const auto begin_itr = lower_bound( disconts, arg_start, [] (const res_arr_idx_pair &p, const res_arrow &a) { return p.first < a; } );
	const auto end_itr   = upper_bound( disconts, arg_stop,  [] (const res_arrow &a, const res_arr_idx_pair &p) { return a < p.first; } );
	return irange(
		numeric_cast<size_t>( distance( common::cbegin( disconts ), begin_itr ) ),
		numeric_cast<size_t>( distance( common::cbegin( disconts ), end_itr   ) )
	);
}

/// \brief Get the size (ie number of discontiguous hits)
size_t discont_hits_index_by_start::size() const {
	return disconts.size();
}

/// \brief Return a reference to the calc_hit associated with an index
///
/// The only source of sensible arguments to this method is get_index_indices_of_disconts_in_range()
const calc_hit & discont_hits_index_by_start::get_discont_hit_of_index_index(const size_t &arg_index_index ///< The index (in this discont_hits_index_by_start) of the calc_hit to return
                                                                             ) const {
	return the_hits.get()[ disconts[ arg_index_index ].second ];
}

/// \brief Generate a string describing the specified discont_hits_index_by_start
///
/// \relates discont_hits_index_by_start
string cath::rslv::detail::to_string(const discont_hits_index_by_start &arg_dhibs ///< The discont_hits_index_by_start to describe
                                     ) {
	return 
		"discont_hits_index_by_start[ "
		+ join(
			irange( 0_z, arg_dhibs.size() )
				| transformed( [&] (const size_t &x) { return to_string( arg_dhibs.get_discont_hit_of_index_index( x ) ); } ),
			", "
		)
		+ " ]";
}