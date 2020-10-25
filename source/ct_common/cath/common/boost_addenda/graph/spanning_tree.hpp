/// \file
/// \brief The spanning_tree class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_GRAPH_SPANNING_TREE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_GRAPH_SPANNING_TREE_HPP

#include <boost/filesystem/path.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/boost_addenda/range/stable_sort_proj.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/common/type_traits.hpp"

#include <tuple>
#include <vector>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Type alias for the Boost Graph edge_descriptor of a graph type
			template <typename G>
			using edge_descriptor_t = typename boost::graph_traits<G>::edge_descriptor;

			/// \brief Type alias for a tuple of the values of an edge range and a weight range
			template <typename EdgeRng,
			          typename WeightRng>
			using edge_wght_tpl_t = decltype( std::tuple_cat(
				std::declval<             range_value_t< EdgeRng   >   >(),
				std::declval< std::tuple< range_value_t< WeightRng > > >()
			) );

			/// \brief Type alias for a vector of edge_wght_tpl_t<>
			template <typename EdgeRng,
			          typename WeightRng>
			using edge_wght_tpl_vec_t = std::vector<edge_wght_tpl_t<EdgeRng, WeightRng>>;

			/// \brief Represent whether a spanning tree should be a max or min spanning tree 
			enum class spanning_tree_dirn : bool {
				MAX, ///< A max spanning tree
				MIN  ///< A min spanning tree
			};

			/// \brief Implementation function to calculate a spanning tree from the specified range of edges and range of weights
			///        for the specified number of items
			template <spanning_tree_dirn Dirn,
			          typename EdgeRng,
			          typename WeightRng>
			edge_wght_tpl_vec_t<EdgeRng, WeightRng> calc_spanning_tree_impl(const EdgeRng   &prm_edges,    ///< The range of edges
			                                                                const WeightRng &prm_weights,  ///< The range of weights
			                                                                const size_t    &prm_num_items ///< The number of items to be spanned
			                                                                ) {
				// If there are zero/one items, just return empty
				if ( prm_num_items <= 1 ) {
					return {};
				}

				// Prepare some type aliases that are useful for this
				using graph     = boost::adjacency_list<
					boost::vecS,
					boost::vecS,
					boost::undirectedS,
					boost::no_property,
					boost::property< boost::edge_weight_t, range_value_t<WeightRng> >
				>;
				using edge_desc = edge_descriptor_t<graph>;

				const auto transformed_weights = prm_weights | boost::adaptors::transformed( [] (const range_value_t< WeightRng> &x) {
					return ( Dirn == spanning_tree_dirn::MAX ) ? -x : x;
				} );

				// Construct a graph from the edges and weights
				const graph my_graph(
					common::cbegin( prm_edges           ),
					common::cend  ( prm_edges           ),
					common::cbegin( transformed_weights ),
					prm_num_items
				);

				// Call kruskal_minimum_spanning_tree() to construct the spanning tree
				std::vector<edge_desc> spanning_tree;
				spanning_tree.reserve( prm_num_items );
				kruskal_minimum_spanning_tree( my_graph, back_inserter( spanning_tree ) );

				// If the number of results edges isn't one less than the number of items, then it was not
				// possible to find a single tree to span all items so throw an exception
				if ( spanning_tree.size() + 1 != prm_num_items ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unable to find tree to span all items"));
				}

				// Build results and then return the result of stably-sorting them
				return stable_sort_proj_copy(
					transform_build<edge_wght_tpl_vec_t<EdgeRng, WeightRng>>(
						spanning_tree,
						[&] (const edge_desc &spanning_tree_edge) {
							const auto &the_weight = get( boost::edge_weight, my_graph, spanning_tree_edge );
							return std::make_tuple(
								source( spanning_tree_edge, my_graph ),
								target( spanning_tree_edge, my_graph ),
								( Dirn == spanning_tree_dirn::MAX ) ? -the_weight : the_weight
							);
						}
					),
					std::less<>{},
					[&] (const edge_wght_tpl_t<EdgeRng, WeightRng> &x) {
						constexpr size_t  x_tpl_size = std::tuple_size< common::remove_cvref_t< decltype( x ) > >::value;
						const     auto   &last_value = std::get< x_tpl_size - 1>( x );
						return std::make_tuple(
							( Dirn == spanning_tree_dirn::MAX ) ? -last_value : last_value,
							std::get< 0 >( x ),
							std::get< 1 >( x )
						);
					}
				);
			}

		} // namespace detail

		/// \brief Calculate a min spanning tree from the specified range of edges and range of weights
		///        for the specified number of items
		template <typename EdgeRng,
		          typename WeightRng>
		detail::edge_wght_tpl_vec_t<EdgeRng, WeightRng> calc_min_spanning_tree(const EdgeRng   &prm_edges,    ///< The range of edges
		                                                                       const WeightRng &prm_weights,  ///< The range of weights
		                                                                       const size_t    &prm_num_items ///< The number of items to be spanned
		                                                                       ) {
			return detail::calc_spanning_tree_impl<detail::spanning_tree_dirn::MIN>(
				prm_edges,
				prm_weights,
				prm_num_items
			);
		}

		/// \brief Calculate a max spanning tree from the specified range of edges and range of weights
		///        for the specified number of items
		template <typename EdgeRng,
		          typename WeightRng>
		detail::edge_wght_tpl_vec_t<EdgeRng, WeightRng> calc_max_spanning_tree(const EdgeRng   &prm_edges,    ///< The range of edges
		                                                                       const WeightRng &prm_weights,  ///< The range of weights
		                                                                       const size_t    &prm_num_items ///< The number of items to be spanned
		                                                                       ) {
			return detail::calc_spanning_tree_impl<detail::spanning_tree_dirn::MAX>(
				prm_edges,
				prm_weights,
				prm_num_items
			);
		}

		size_size_pair_vec make_simple_unweighted_spanning_tree(const size_t &);
		size_size_doub_tpl_vec calc_max_spanning_tree(const size_size_doub_tpl_vec &,
		                                              const size_t &);
		size_size_doub_tpl_vec calc_min_spanning_tree(const size_size_doub_tpl_vec &,
		                                              const size_t &);
		size_size_pair_vec get_edges_of_spanning_tree(const size_size_doub_tpl_vec &);
		size_size_doub_tpl_vec order_spanning_tree_from_start(const size_size_doub_tpl_vec &,
		                                                      const size_t &);
		std::string make_graphviz_string_of_spanning_tree(const size_size_doub_tpl_vec &);
		void write_graphviz_string_of_spanning_tree(const size_size_doub_tpl_vec &,
		                                            const boost::filesystem::path &);

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_GRAPH_SPANNING_TREE_HPP
