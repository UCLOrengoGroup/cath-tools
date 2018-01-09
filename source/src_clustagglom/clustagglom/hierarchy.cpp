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
#include "clustagglom/links.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/back.hpp"
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

/// \brief Get the indices of the representatives in the clusters of the specified hierarchy
///
/// \pre arg_hierarchy must be exactly two layers deep (ie corresponding to one level of clustering)
///
/// \relates hierarchy
size_vec cath::clust::get_rep_indices(const hierarchy &arg_hierarchy ///< The hierarchy to query
                                      ) {
	if ( arg_hierarchy.size() != 2 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to get reps from clustering hierarchy that doesn't have exactly one level of clustering"));
	}

	return transform_build<size_vec>(
		front( front( arg_hierarchy ) ),
		[&] (const hierarchy_value &x) {
			return front( back( arg_hierarchy )[ x.get_index() ] ).get_index();
		}
	);
}

/// \brief Get the indices of the groups of items in the clusters of the specified hierarchy
///
/// \pre arg_hierarchy must be exactly two layers deep (ie corresponding to one level of clustering)
///
/// \relates hierarchy
size_vec_vec cath::clust::get_index_groups(const hierarchy &arg_hierarchy ///< The hierarchy to query
                                           ) {
	if ( arg_hierarchy.size() != 2 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to get reps from clustering hierarchy that doesn't have exactly one level of clustering"));
	}

	return transform_build<size_vec_vec>(
		front( front( arg_hierarchy ) ),
		[&] (const hierarchy_value &x) {
			return ( x.get_type() == hierarchy_ref::ENTRY )
				? size_vec{ x.get_index() }
				: transform_build<size_vec>(
					back( arg_hierarchy )[ x.get_index() ],
					[] (const hierarchy_value &y) {
						return y.get_index();
					}
				);
		}
	);
}

/// \brief Get spanning trees for the groups of items in the clusters of the specified hierarchy
///
/// \pre arg_hierarchy must be exactly two layers deep (ie corresponding to one level of clustering)
///
/// \relates hierarchy
size_size_pair_vec_vec cath::clust::get_spanning_trees(const hierarchy &arg_hierarchy, ///< The hierarchy to query
                                                       const links     &arg_links      ///< The links between the items
                                                       ) {
	return transform_build<size_size_pair_vec_vec>(
		get_index_groups( arg_hierarchy ),
		[&] (const size_vec &index_group) {
			return get_spanning_tree_of_subset(
				arg_links,
				size_set{ common::cbegin( index_group ), common::cend( index_group ) }
			);
		}
	);
}

/// \brief Write the specified spanning trees to the specified stream using the specified names
ostream & cath::clust::write_spanning_trees(ostream                      &arg_os,       ///< The ostream to which the spanning trees should be written
                                            const size_size_pair_vec_vec &arg_trees,    ///< The spanning trees to write
                                            const id_of_str_bidirnl      &arg_name_ider ///< The holding the IDs of the entries in the hierarchy
                                            ) {
	for (const size_size_pair_vec &tree : arg_trees) {
		for (const size_size_pair &edge : tree) {
			arg_os
				<< arg_name_ider.get_name_of_id( edge.first  )
				<< " "
				<< arg_name_ider.get_name_of_id( edge.second )
				<< "\n";
		}
	}
	return arg_os;
}

/// \brief Calculate spanning trees for the clusters in the specified hierarchy and write them to the specified ostream
///
/// \pre arg_hierarchy must be exactly two layers deep (ie corresponding to one level of clustering)
///
/// \relates hierarchy
ostream & cath::clust::write_spanning_trees(ostream                 &arg_os,        ///< The ostream to which the spanning trees should be written
                                            const hierarchy         &arg_hierarchy, ///< The hierarchy to query
                                            const id_of_str_bidirnl &arg_name_ider, ///< The holding the IDs of the entries in the hierarchy
                                            const links             &arg_links      ///< The links between the items
                                            ) {
	return write_spanning_trees(
		arg_os,
		get_spanning_trees( arg_hierarchy, arg_links ),
		arg_name_ider
	);
}

/// \brief Calculate spanning trees for the clusters in the specified hierarchy and write them to the specified file
///
/// \relates hierarchy
void cath::clust::write_spanning_trees(const path              &arg_output_file, ///< The file to which the reps should be written
                                       const hierarchy         &arg_hierarchy,   ///< The hierarchy to query
                                       const id_of_str_bidirnl &arg_name_ider,   ///< The holding the IDs of the entries in the hierarchy
                                       const links             &arg_links        ///< The links between the items
                                       ) {
	ofstream span_ostream;
	open_ofstream( span_ostream, arg_output_file );
	write_spanning_trees( span_ostream, arg_hierarchy, arg_name_ider, arg_links );
	span_ostream.close();
}

/// \brief Write the reps for the clusters in the specified hierarchy to the specified ostream
///
/// \relates hierarchy
std::ostream & cath::clust::write_reps(ostream                 &arg_os,        ///< The ostream to which the hierarchy should be written
                                       const hierarchy         &arg_hierarchy, ///< The hierarchy to write
                                       const id_of_str_bidirnl &arg_name_ider  ///< The holding the IDs of the entries in the hierarchy
                                       ) {
	// const size_vec rep_indices = get_rep_indices( arg_hierarchy );
	for (const size_t &rep_index : get_rep_indices( arg_hierarchy ) ) {
		arg_os << arg_name_ider.get_name_of_id( rep_index ) << "\n";
	}
	return arg_os;
}

/// \brief Write the reps for the clusters in the specified hierarchy to the specified file
///
/// \relates hierarchy
void cath::clust::write_reps(const path              &arg_output_file, ///< The file to which the reps should be written
                             const hierarchy         &arg_hierarchy,   ///< The hierarchy to write
                             const id_of_str_bidirnl &arg_name_ider    ///< The holding the IDs of the entries in the hierarchy
                             ) {
	ofstream reps_ostream;
	open_ofstream( reps_ostream, arg_output_file );
	write_reps( reps_ostream, arg_hierarchy, arg_name_ider );
	reps_ostream.close();
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
