/// \file
/// \brief The spanning_tree class definitions

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

#include "spanning_tree.hpp"

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>

#include "common/size_t_literal.hpp"

using namespace cath;
using namespace cath::common;

using boost::adaptors::transformed;
using boost::irange;
using std::make_pair;

/// \brief Get a simple ( 0 <-> 1 <-> 2 <-> ... <-> (arg_num_items -1) spanning tree)
size_size_pair_vec cath::common::make_simple_unweighted_spanning_tree(const size_t &arg_num_items ///< The number of items to span
                                                                      ) {
	// If zero/one items then return empty
	if ( arg_num_items <= 1 ) {
		return {};
	}

	// Otherwise return a spanning tree between adjacent indices
	return transform_build<size_size_pair_vec>(
		irange( 1_z, arg_num_items ),
		[] (const size_t &x) { return make_pair( x - 1, x ); }
	);
}

/// \brief Calculate a max-spanning-tree over the specified edges_and_scores and number of items
size_size_doub_tpl_vec cath::common::calc_max_spanning_tree(const size_size_doub_tpl_vec &arg_edges,    ///< The weighted edges from which a spanning tree should be extracted
                                                            const size_t                 &arg_num_items ///< The number of items to span
                                                            ) {
	return calc_max_spanning_tree(
		arg_edges | transformed( [] (const size_size_doub_tpl &x) { return make_pair( get<0>( x ), get<1>( x ) ); } ),
		arg_edges | transformed( [] (const size_size_doub_tpl &x) { return get<2>( x );                           } ),
		arg_num_items
	);
}

/// \brief Calculate a min-spanning-tree over the specified edges_and_scores and number of items
size_size_doub_tpl_vec cath::common::calc_min_spanning_tree(const size_size_doub_tpl_vec &arg_edges,    ///< The weighted edges from which a spanning tree should be extracted
                                                            const size_t                 &arg_num_items ///< The number of items to span
                                                            ) {
	return calc_min_spanning_tree(
		arg_edges | transformed( [] (const size_size_doub_tpl &x) { return make_pair( get<0>( x ), get<1>( x ) ); } ),
		arg_edges | transformed( [] (const size_size_doub_tpl &x) { return get<2>( x );                           } ),
		arg_num_items
	);
}

/// \brief Get the edges of the specified spanning tree (ie strip off the weights)
size_size_pair_vec cath::common::get_edges_of_spanning_tree(const size_size_doub_tpl_vec &arg_spanning_tree ///< The weighted spanning-tree from which the weights should be stripped
                                                            ) {
	return transform_build<size_size_pair_vec>(
		arg_spanning_tree,
		[] (const size_size_doub_tpl &x) {
			return make_pair( get<0>( x ), get<1>( x ) );
		}
	);
}

