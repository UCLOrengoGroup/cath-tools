/// \file
/// \brief The alignment action header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGNMENT_ACTION_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGNMENT_ACTION_HPP

#include <functional>

#include "cath/alignment/align_type_aliases.hpp"
#include "cath/alignment/alignment.hpp"
#include "cath/alignment/aln_glue_style.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/graph/spanning_tree.hpp"
#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/common/boost_addenda/range/max_proj_element.hpp"
#include "cath/common/cpp17/as_const.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_list.hpp"

namespace cath { class protein_list; }
namespace cath { namespace align { class alignment; } }
namespace cath { namespace file { class pdb_list; } }

namespace cath {
	namespace align {
		namespace detail {
			using aln_ent_ind_tup      = std::tuple<const alignment &, const size_t &, const size_t &>;

			using aln_ent_ind_tup_pair = std::pair<aln_ent_ind_tup, aln_ent_ind_tup>;

			enum class glued_row_type : char {
				FROM_A,
				FROM_B,
				FROM_BOTH
			};

			void append_glued_row(alignment &,
			                      const aln_ent_ind_tup_pair &,
			                      const glued_row_type &);
		} // namespace detail

		alignment glue_two_alignments(const alignment &,
		                              const size_t &,
		                              const alignment &,
		                              const size_t & = 0);

		alignment build_alignment_from_parts(const size_size_alignment_tuple_vec &,
		                                     const protein_list &,
		                                     const aln_glue_style &);

		/// \brief Build an alignment between the specified proteins with specified similarity scores with the specified approach
		///        using the specified function to get pairwise alignments
		template <typename Fn>
		std::pair<alignment, size_size_doub_tpl_vec> build_alignment(const protein_list            &prm_proteins, ///< The proteins for which the alignment is to be built
		                                                             const size_size_doub_tpl_vec  &prm_scores,   ///< The scores between the proteins in the alignment
		                                                             const aln_glue_style          &prm_strategy, ///< The approach that should be used for glueing alignments together
		                                                             Fn                           &&prm_fn        ///< A callback function for getting the alignment between the two proteins of the specified indices
		                                                             ) {
			// If there's only zero/one protein(s), just return an appropriate simple alignment
			if ( prm_proteins.empty() ) {
				return { alignment{ 0 }, size_size_doub_tpl_vec{} };
			}
			if ( prm_proteins.size() == 1 ) {
				return {
					make_single_alignment( common::front( prm_proteins ).get_length() ),
					size_size_doub_tpl_vec{}
				};
			}

			// Get a spanning tree
			const auto spanning_tree = [&] {
				auto lcl_span_tree = common::calc_max_spanning_tree( prm_scores, prm_proteins.size() );

				// If aln_glue_style::INCREMENTALLY_WITH_PAIR_REFINING, order spanning tree
				if ( prm_strategy == aln_glue_style::INCREMENTALLY_WITH_PAIR_REFINING ) {
					// At present, just start incremental building from the highest-scoring entry in the spanning tree
					//
					// It might be better if this was chosen by a process like: repeatedly
					// find the worst score and then restrict to larger side of the remaining
					// tree until there's only one link left
					const size_t START_INDEX = boost::numeric_cast<size_t>( std::distance(
						common::cbegin( lcl_span_tree ),
						common::max_proj_element(
							common::as_const( lcl_span_tree ),
							std::less<>{},
							[] (const size_size_doub_tpl &x) { return get<2>( x ); }
						)
					) );
					lcl_span_tree = common::order_spanning_tree_from_start( lcl_span_tree, START_INDEX );
				}
				return lcl_span_tree;
			} ();

			// Get alignments
			const auto alignments_tree = common::transform_build<size_size_alignment_tuple_vec>(
				spanning_tree,
				[&] (const size_size_doub_tpl &spanning_tree_edge) {
					const size_t &index_a = std::get<0>( spanning_tree_edge );
					const size_t &index_b = std::get<1>( spanning_tree_edge );
					return std::make_tuple(
						index_a,
						index_b,
						::std::invoke(
							prm_fn,
							index_a,
							index_b
						)
					);
				}
			);

			// Make with build_alignment_from_parts
			return make_pair(
				build_alignment_from_parts(
					alignments_tree,
					prm_proteins,
					prm_strategy
				),
				spanning_tree
			);
		}

	} // namespace align
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGNMENT_ACTION_HPP
