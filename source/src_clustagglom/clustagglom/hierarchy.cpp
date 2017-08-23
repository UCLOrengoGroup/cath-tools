/// \file
/// \brief The hierarchy class definitions

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

#include "hierarchy.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/reverse.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "clustagglom/hierarchy/hierarchy_fn.hpp"
#include "clustagglom/hierarchy/hierarchy_value.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/container/id_of_str_bidirnl.hpp"
#include "common/file/open_fstream.hpp"
#include "common/size_t_literal.hpp"

#include <fstream>
#include <string>

using namespace cath;
using namespace cath::clust;
using namespace cath::common;

using boost::adaptors::reversed;
using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::filesystem::path;
using boost::format;
using boost::range::reverse;
using boost::range::sort;
using std::ofstream;
using std::ostream;
using std::string;

/// \brief Get the index of the first entry in the specified hierarchy
///        from the specified value at the specified depth
///
/// If the specified value is an entry type, this just returns that value's index
/// otherwise it calls itself one layer deeper.
///
/// This is useful in sort_hierarchy()
///
/// \relates hierarchy
item_idx cath::clust::detail::first_index(const hierarchy       &arg_hierarchy,  ///< The hierarchy to search
                                          const hierarchy_value &arg_value,      ///< The value from which to start looking
                                          const size_t          &arg_value_depth ///< The depth of the value
                                          ) {
	return ( arg_value.get_type() == hierarchy_ref::CLUSTER )
		? first_index(
			arg_hierarchy,
			front( arg_hierarchy[ arg_value_depth + 1 ][ arg_value.get_index() ] ),
			arg_value_depth + 1
		)
		: arg_value.get_index();
}

/// \brief Generate a string describing the specified hierarchy
///
/// \relates hierarchy
string cath::clust::to_string(const hierarchy &arg_hierarchy ///< The hierarchy to describe
                              ) {
	using std::to_string;
	return
		  "hierarchy["
		+ join(
			indices( arg_hierarchy.size() )
				| transformed( [&] (const size_t &x) {
					return
						  "\n\tLAYER "
						+ to_string( x )
						+ ":\n\t\t"
						+ to_string( arg_hierarchy[ x ], "\n\t\t" )
						+ "\n";
				} ),
			""
		)
		+ "]";
}

/// \brief Insert a description of the specified hierarchy into the specified ostream
///
/// \relates hierarchy
std::ostream & cath::clust::operator<<(ostream         &arg_os,       ///< The ostream into which the description should be inserted
                                       const hierarchy &arg_hierarchy ///< The hierarchy to describe
                                       ) {
	arg_os << to_string( arg_hierarchy );
	return arg_os;
}

/// \brief Sort the specified hierarchy according to the specified sorting indices
///
/// \relates hierarchy
void cath::clust::sort_hierarchy(hierarchy      &arg_hierarchy,      ///< The hierarchy to sort
                                 const size_vec &arg_sorting_indices ///< Values corresponding to each of the entries such that one value is less than another if the corresponding entry should be sorted before the other
                                 ) {
	for (const size_t &depth_idx : indices( arg_hierarchy.size() ) | reversed ) {
		hierarchy_layer &layer = arg_hierarchy[ depth_idx ];
		for (hierarchy_group &group : layer) {
			sort(
				group,
				[&] (const hierarchy_value &x, const hierarchy_value &y) {
					return (
						arg_sorting_indices[ detail::first_index( arg_hierarchy, x, depth_idx ) ]
						<
						arg_sorting_indices[ detail::first_index( arg_hierarchy, y, depth_idx ) ]
					);
				}
			);
		}
	}
}

/// \brief Copy the specified hierarchy and then sort they copy according to the specified sorting indices and return it
///
/// \relates hierarchy
hierarchy cath::clust::sort_hierarchy_copy(hierarchy       arg_hierarchy,      ///< The hierarchy to sort
                                           const size_vec &arg_sorting_indices ///< Values corresponding to each of the entries such that one value is less than another if the corresponding entry should be sorted before the other
                                           ) {
	sort_hierarchy( arg_hierarchy, arg_sorting_indices );
	return arg_hierarchy;
}

/// \brief Write the clusters of the specified hierarchy to the specified ostream
///        using the IDs in the specified IDer
///
/// \pre The IDer must correspond to the hierarchy
///      (and therefore have an ID for each entry in the hierarchy)
///
/// \relates hierarchy
ostream & cath::clust::write_cluster(ostream                 &arg_os,        ///< The ostream to which the hierarchy should be written
                                     const hierarchy         &arg_hierarchy, ///< The hierarchy to write
                                     const id_of_str_bidirnl &arg_name_ider  ///< The holding the IDs of the entries in the hierarchy
                                     ) {
	using std::to_string;
	detail::depth_first_traverse_hierachy(
		arg_hierarchy,
		[&] (const item_vec &counters,
		     const item_idx &entry_index
		     ) {
			arg_os
				<< arg_name_ider.get_name_of_id( entry_index )
				<< join(
					counters
						| transformed( [] (const item_idx &x) { return "\t" + to_string( x );} ),
					""
				)
				<< "\n";
		}
	);
	return arg_os;
}

/// \brief Write the clusters of the specified hierarchy to the specified file
///        using the IDs in the specified IDer
///
/// \pre The IDer must correspond to the hierarchy
///      (and therefore have an ID for each entry in the hierarchy)
///
/// \relates hierarchy
void cath::clust::write_cluster(const path              &arg_output_file, ///< The file to which the hierarchy should be written
                                const hierarchy         &arg_hierarchy,   ///< The hierarchy to write
                                const id_of_str_bidirnl &arg_name_ider    ///< The holding the IDs of the entries in the hierarchy
                                ) {
	ofstream clust_ostream;
	open_ofstream( clust_ostream, arg_output_file );
	write_cluster( clust_ostream, arg_hierarchy, arg_name_ider );
	clust_ostream.close();
}

/// \brief Make a hierarchy from the specified layers in reversed order (ie deepest first)
///
/// \pre `reverse(arg_layers)` must still make a valid hierarchy
///
/// \relates hierarchy
hierarchy cath::clust::make_hierarchy_from_reversed_layers(hierarchy_layer_vec arg_layers ///< The reversed layers from which to make the hierarchy (ie deepest first)
                                                           ) {
	reverse( arg_layers );
	return hierarchy{ move( arg_layers ) };
}
