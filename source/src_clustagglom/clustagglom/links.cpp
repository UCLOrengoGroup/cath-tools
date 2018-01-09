/// \file
/// \brief The links class definitions

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

#include "links.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/stable_sort.hpp>
#include <boost/range/combine.hpp>

#include "clustagglom/link.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/boost_addenda/graph/spanning_tree.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/container/id_of_str_bidirnl.hpp"
#include "common/file/open_fstream.hpp"
#include "common/size_t_literal.hpp"

#include <fstream>
#include <tuple>

using namespace cath;
using namespace cath::common;
using namespace cath::clust;

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::filesystem::path;
using boost::range::combine;
using boost::range::for_each;
using boost::range::sort;
using boost::range::stable_sort;
using std::get;
using std::ofstream;
using std::ostream;
using std::string;
using std::tie;

/// \brief Make links from the specified raw links data
///
/// \relates links
links cath::clust::make_links(const item_item_strength_tpl_vec &arg_raw_links ///< The raw links data from which the links should be built
                              ) {
	links result;
	for_each(
		arg_raw_links,
		[&] (const item_item_strength_tpl &x) { add_link_symmetrically( result, x ); }
	);
	return result;
}

/// \brief Get a spanning tree for the specified subset of items in the specified links
///
/// \pre The items in arg_index_group must be spanned by arg_links
///
/// \relates links
size_size_pair_vec cath::clust::get_spanning_tree_of_subset(const links    &arg_links,      ///< The links from which the spanning tree should be formed
                                                            const size_set &arg_index_group ///< The items over which the spanning tree should be formed
                                                            ) {
	size_size_doub_tpl_vec relevant_links;
	relevant_links.reserve( arg_index_group.size() );
	for (const size_t &index : arg_index_group) {
		for (const link &x : arg_links[ index ] ) {
			if ( x.node < index && contains( arg_index_group, x.node ) ) {
				relevant_links.emplace_back(
					x.node,
					index,
					x.dissim
				);
			}
		}
	}

	return get_edges_of_spanning_tree( calc_min_spanning_tree( relevant_links, arg_index_group.size() ) );
}

/// \brief Generate a string describing the specified links
///
/// \relates links
string cath::clust::to_string(const links &arg_links ///< The links to describe
                              ) {
	return "links["
		+ join(
			indices( arg_links.size() )
				| transformed( [&] (const size_t &x) {
					return "\n\t" + link_list_string( arg_links[ x ], x );
				} ),
			""
		)
		+ "\n]";
}

/// \brief Write the specified links to the specified ostream using the specified
///        name_ider and sort indices
///
/// The aim of this function is to reorder the links in such a way that TCluster
/// will give the same output as this code (because, unlike this code, TCluster
/// gives slightly different results depending on the order in which it receives links)
///
/// \relates links
void cath::clust::write_ordered_links(ostream                 &arg_output,      ///< The ostream to which the ordered links should be written
                                      const links             &arg_links,       ///< The links to be written
                                      const id_of_str_bidirnl &arg_ider,        ///< The name_ider containing the names for the links
                                      const size_vec          &arg_sort_indices ///< The ranks of the items (ie a 0 should appear in the index corresponding to that of the most preferred item)
                                      ) {
	/// \TODO: Move this to common/type_aliases.hpp
	using size_size_size_size_tpl_vec = std::vector<size_size_size_size_tpl>;

	size_size_size_size_tpl_vec print_links;
	for (const size_t &outer_index : indices( arg_links.size() ) ) {
		const auto &inner_links      = arg_links       [ outer_index ];
		const auto &outer_sort_index = arg_sort_indices[ outer_index ];
		const auto &outer_name       = arg_ider.get_name_of_id( outer_index );
		for (const boost::tuple<const link &, size_t> &x : combine( inner_links, indices( inner_links.size() ) ) ) {
			const auto &inner_link       = x.get<0>();
			const auto &inner_list_index = x.get<1>();
			const auto &inner_index      = inner_link.node;
			if ( outer_name < arg_ider.get_name_of_id( inner_index ) ) {
				const auto &inner_sort_index = arg_sort_indices[ inner_index ];

				print_links.emplace_back(
					outer_index,
					inner_list_index,
					std::max( outer_sort_index, inner_sort_index ),
					std::min( outer_sort_index, inner_sort_index )
				);
			}
		}
	}

	sort(
		print_links,
		[] (const size_size_size_size_tpl &x, const size_size_size_size_tpl &y) {
			return (
				tie( get<3>( x ), get<2>( x ) )
				<
				tie( get<3>( y ), get<2>( y ) )
			);
		}
	);

	for (const auto &print_link : print_links) {
		const auto &the_link = arg_links[ get<0>( print_link ) ][ get<1>( print_link ) ];
		arg_output
			<< arg_ider.get_name_of_id( get<0>( print_link ) )
			<< " "
			<< arg_ider.get_name_of_id( the_link.node )
			<< "\t"
			<< std::setprecision( 30 ) << ( -the_link.dissim )
			<< "\t100\n";
	}
}

/// \brief Write the specified links to the specified file using the specified
///        name_ider and sort indices
///
/// \relates links
void cath::clust::write_ordered_links(const path              &arg_output,      ///< The file to which the ordered links should be written
                                      const links             &arg_links,       ///< The links to be written
                                      const id_of_str_bidirnl &arg_ider,        ///< The name_ider containing the names for the links
                                      const size_vec          &arg_sort_indices ///< The ranks of the items (ie a 0 should appear in the index corresponding to that of the most preferred item)
                                      ) {
	ofstream links_ostream;
	open_ofstream( links_ostream, arg_output );
	write_ordered_links( links_ostream, arg_links, arg_ider, arg_sort_indices );
	links_ostream.close();
}
