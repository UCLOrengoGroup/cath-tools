/// \file
/// \brief The chimera_viewer class definitions

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

#include "chimera_viewer.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "alignment/alignment.hpp"
#include "chopping/region/region.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/batch/batch_functions.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "display/viewer/pymol/pymol_tools.hpp"
#include "display_colour/display_colour.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "superposition/io/superposition_io.hpp"
#include "superposition/superposition_context.hpp"

#include <algorithm>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::detail;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sup;

using boost::adaptors::transformed;
using boost::adjacency_list;
using boost::algorithm::join;
using boost::algorithm::replace_all_copy;
using boost::edge_weight_t;
using boost::graph_traits;
using boost::irange;
using boost::lexical_cast;
using boost::no_property;
using boost::property;
using boost::string_ref;
using boost::undirectedS;
using boost::vecS;
using std::make_pair;
using std::map;
using std::min;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::string;
using std::vector;

constexpr size_t chimera_viewer::RESIDUE_BATCH_SIZE;

/// \brief TODOCUMENT
///
/// \relates chimera_viewer
///
/// \todo Refactor out any similarities between write_chimera_pair_alignments() and write_chimera_global_alignment()
///
/// \todo Separate out the code that identifies the residue links in the alignment from the code that writes
///       this information in Chimera-ese (so that the first bit of code can be reused with other viewers)
void cath::detail::write_chimera_pair_alignments(ostream                     &prm_os,                   ///< TODOCUMENT
                                                 const superposition_context &prm_superposition_context ///< TODOCUMENT
                                                 ) {
	if ( ! prm_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot write PyMOL pair alignments for superposition_context with no alignment"));
	}
	const alignment     &the_alignment     = prm_superposition_context.get_alignment();
// 	const superposition &the_superposition = prm_superposition_context.get_superposition();
	const pdb_list       pdbs              = get_restricted_pdbs( prm_superposition_context );
	const str_vec        names             = clean_names_for_viewer( prm_superposition_context );
	
	// Grab some basic details
	const alignment::size_type num_entries = min( the_alignment.num_entries(), names.size() );
	const alignment::size_type aln_length  = the_alignment.length();
	const residue_id_vec_vec   residue_ids = get_backbone_complete_residue_ids_of_first_chains( pdbs );

	for (const size_t &entry_ctr_a : indices( num_entries ) ) {
		const string         &name_a        = names      [ entry_ctr_a ];
		const residue_id_vec &residue_ids_a = residue_ids[ entry_ctr_a ];

		for (const size_t &entry_ctr_b : irange( entry_ctr_a + 1, num_entries ) ) {
			const string         &name_b        = names      [ entry_ctr_b ];
			const residue_id_vec &residue_ids_b = residue_ids[ entry_ctr_b ];

			bool added_pair_distances(false);
			for (const alignment::size_type &aln_posn_ctr : indices( aln_length ) ) {
				const aln_posn_opt position_a = the_alignment.position_of_entry_of_index( entry_ctr_a, aln_posn_ctr );
				const aln_posn_opt position_b = the_alignment.position_of_entry_of_index( entry_ctr_b, aln_posn_ctr );
				if ( position_a && position_b) {
					added_pair_distances = true;
					if ( *position_a >= residue_ids_a.size() ) {
						BOOST_THROW_EXCEPTION(invalid_argument_exception(
							"Whilst adding alignment extras in chimera_viewer, residue index "
							+ lexical_cast<string>( *position_a )
							+ " is out of range "
							+ lexical_cast<string>(residue_ids_a.size())
						));
					}
					if ( *position_b >= residue_ids_b.size()) {
						BOOST_THROW_EXCEPTION(invalid_argument_exception(
							"Whilst adding alignment extras in chimera_viewer, residue index "
							+ lexical_cast<string>( *position_b )
							+ " is out of range "
							+ lexical_cast<string>(residue_ids_b.size())
						));
					}
					const residue_id &residue_id_a = residue_ids_a[ *position_a ];
					const residue_id &residue_id_b = residue_ids_b[ *position_b ];

					prm_os << "distance ";
					prm_os << name_a;
					prm_os << "_";
					prm_os << name_b;
					prm_os << "_alignment, /";
					prm_os << name_a;
					prm_os << "///";
					prm_os << chimera_viewer::parse_residue_id_for_chimera( residue_id_a );
					prm_os << "/CA/, /";
					prm_os << name_b;
					prm_os << "///";
//					prm_os << residue_id_b;
					prm_os << chimera_viewer::parse_residue_id_for_chimera( residue_id_b );
					prm_os << "/CA/\n";
				}
			}

			if (added_pair_distances) {
				prm_os << "disable ";
				prm_os << name_a;
				prm_os << "_";
				prm_os << name_b;
				prm_os << "_alignment\n";
			}
		}
	}
	prm_os << "hide labels\n";
	prm_os << "set dash_gap,    0.0\n";
	prm_os << "set dash_color,  black\n";
	prm_os << "set dash_radius, 0.05\n";

}

/// \brief TODOCUMENT
///
/// \relates chimera_viewer
///
/// \todo Refactor out any similarities between write_chimera_pair_alignments() and write_chimera_global_alignment()
///
/// \todo Separate out the code that identifies the residue links in each alignment from the code that writes
///       this information in Chimera-ese (so that the first bit of code can be reused with other viewers)
void cath::detail::write_chimera_global_alignment(ostream                     &prm_os,                   ///< TODOCUMENT
                                                  const superposition_context &prm_superposition_context ///< TODOCUMENT
                                                  ) {
	if ( ! prm_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot write PyMOL global alignment for superposition_context with no alignment"));
	}
	const alignment     &the_alignment     = prm_superposition_context.get_alignment();
	const superposition &the_superposition = prm_superposition_context.get_superposition();
	const pdb_list       pdbs              = get_restricted_pdbs( prm_superposition_context );
	const str_vec        names             = clean_names_for_viewer( prm_superposition_context );
	
	// Grab some basic details
	const alignment::size_type num_entries = min( the_alignment.num_entries(), names.size() );
	const alignment::size_type aln_length  = the_alignment.length();
	const residue_id_vec_vec   residue_ids = get_backbone_complete_residue_ids_of_first_chains( pdbs );

	/// ???
	bool added_distances(false);
	for (const alignment::size_type &aln_index : indices( aln_length ) ) {
		// Prepare some type aliases that are useful for this
		using Graph     = adjacency_list<vecS, vecS, undirectedS, no_property, property <edge_weight_t, double>>;
		using edge_desc = graph_traits<Graph>::edge_descriptor;

		// Prepare edges and distances to be loaded into the graph
		size_size_pair_vec edges;
		doub_vec           distances;

		// ????
		const size_t num_present_posns = num_present_positions_of_index( the_alignment, aln_index );
		for (const size_t &entry_a : indices( num_entries ) ) {
			for (const size_t &entry_b : irange( entry_a + 1, num_entries ) ) {

				const aln_posn_opt posn_a = the_alignment.position_of_entry_of_index( entry_a, aln_index );
				const aln_posn_opt posn_b = the_alignment.position_of_entry_of_index( entry_b, aln_index );

				if ( posn_a && posn_b ) {
					added_distances = true;

//					cerr << "(1) At alignment index\t" << aln_index << ",\tabout to get CA coord from\t" << entry_a << "\t(\t" << names[ entry_a ] << "\t), index:\t" << *posn_a << endl;
//					cerr << "(2) At alignment index\t" << aln_index << ",\tabout to get CA coord from\t" << entry_b << "\t(\t" << names[ entry_b ] << "\t), index:\t" << *posn_b << endl;
					const coord         ca_a        = get_residue_ca_coord_of_backbone_complete_index( pdbs[ entry_a ], *posn_a );
					const coord         ca_b        = get_residue_ca_coord_of_backbone_complete_index( pdbs[ entry_b ], *posn_b );
					const double        distance    = superposed_distance( the_superposition, entry_a, ca_a, entry_b, ca_b );

					edges.push_back    ( make_pair( entry_a, entry_b ) );
					distances.push_back( distance                      );
				}
			}
		}

		// Construct a graph from the edges and weights
		const Graph my_graph(
			common::cbegin( edges     ),
			common::cend  ( edges     ),
			common::cbegin( distances ),
			num_present_posns
		);

		// Call kruskal_minimum_spanning_tree() to construct the spanning tree
		vector<edge_desc> spanning_tree;
		kruskal_minimum_spanning_tree( my_graph, back_inserter( spanning_tree ) );

		// ????
		for (const edge_desc &spanning_tree_edge : spanning_tree) {
			const size_t         entry_a     = source( spanning_tree_edge, my_graph );
			const size_t         entry_b     = target( spanning_tree_edge, my_graph );
			const aln_posn_type  res_index_a = get_position_of_entry_of_index( the_alignment, entry_a, aln_index );
			const aln_posn_type  res_index_b = get_position_of_entry_of_index( the_alignment, entry_b, aln_index );
			const string        &name_a      = names    [ entry_a ];
			const string        &name_b      = names    [ entry_b ];
			const residue_id    &res_name_a  = residue_ids[ entry_a ][ res_index_a ];
			const residue_id    &res_name_b  = residue_ids[ entry_b ][ res_index_b ];

			prm_os << "distance ";
			prm_os << "alignment, /";
			prm_os << name_a;
			prm_os << "///";
			prm_os << chimera_viewer::parse_residue_id_for_chimera( res_name_a );
			prm_os << "/CA/, /";
			prm_os << name_b;
			prm_os << "///";
			prm_os << chimera_viewer::parse_residue_id_for_chimera( res_name_b );
			prm_os << "/CA/\n";
		}
	}
	if (added_distances) {
		prm_os << "disable alignment\n";
	}

	enum class coreness : bool {
		NONCORE,
		CORE,
	};
	using str_res_id_vec_map              = map<string, residue_id_vec>;
	using coreness_str_res_id_vec_map_map = map<coreness, str_res_id_vec_map>;

	// ????
	if ( the_alignment.is_scored() ) {
		coreness_str_res_id_vec_map_map core_res_ids_of_entry_name;
		const alignment_residue_scores &the_scores = the_alignment.get_alignment_residue_scores();
		for (const size_t &entry : indices( num_entries ) ) {
			const string &entry_name = names[ entry ];
			for (const alignment::size_type &index : indices( aln_length ) ) {
				if ( has_score( the_scores, entry, index ) ) {
					const coreness is_core = ( get_score( the_scores, entry, index, true, true ) > 0.25 )
					                         ? coreness::CORE
					                         : coreness::NONCORE;
					core_res_ids_of_entry_name[ is_core ][ entry_name ].push_back(
						residue_ids[ entry ][ get_position_of_entry_of_index( the_alignment, entry, index ) ]
					);
				}
			}
		}
		prm_os << join(
			core_res_ids_of_entry_name
				// \TODO Come C++17 and structured bindings, use here
				| transformed( [] (const pair<const coreness, str_res_id_vec_map> &core_data) {
					const coreness           &is_core               = core_data.first;
					const str_res_id_vec_map &res_ids_of_entry_name = core_data.second;

					// \TODO Come C++17 and structured bindings, use here
					str_vec selection_strings;
					for (const pair<const string, residue_id_vec> &entry_name_and_res_ids : res_ids_of_entry_name) {
						const string         &entry_name = entry_name_and_res_ids.first;
						const residue_id_vec &res_ids    = entry_name_and_res_ids.second;

						const size_t num_res_ids     = res_ids.size();
						const size_t num_res_batches = num_batches( num_res_ids, chimera_viewer::RESIDUE_BATCH_SIZE, broken_batch_tol::PERMIT );
						for (const size_t &batch_ctr : indices( num_res_batches ) ) {
							const size_size_pair begin_and_end = batch_begin_and_end(
								num_res_ids,
								chimera_viewer::RESIDUE_BATCH_SIZE,
								batch_ctr,
								broken_batch_tol::PERMIT
							);
							selection_strings.push_back( pymol_tools::pymol_res_seln_str(
								entry_name,
								residue_id_vec{
									std::next( common::cbegin( res_ids ), static_cast<ptrdiff_t>( begin_and_end.first  ) ),
									std::next( common::cbegin( res_ids ), static_cast<ptrdiff_t>( begin_and_end.second ) )
								}
							) );
						}
					}
					return
						  "select "
						+ ( ( is_core == coreness::CORE ) ? "core"s : "noncore"s )
						+ ", ( "
						+ join( selection_strings, " or " )
						+ " )";
				} ),
			"\n"
		) << "\n";
		prm_os << "deselect\n";
	}

	prm_os << "hide labels\n";
	prm_os << "set dash_gap, 0.0\n";
	prm_os << "set dash_color, black\n";
	prm_os << "set dash_radius, 0.05\n";
}


/// \brief TODOCUMENT
string chimera_viewer::do_default_executable() const {
	return "chimera";
}

/// \brief TODOCUMENT
string chimera_viewer::do_default_file_extension() const {
	return ".pml";
}

/// \brief TODOCUMENT
void chimera_viewer::do_write_start(ostream &prm_os ///< TODOCUMENT
                                    ) const {
	prm_os << "feedback disable,all,output\n";
}

/// \brief TODOCUMENT
void chimera_viewer::do_write_load_pdbs(ostream             &prm_os,            ///< TODOCUMENT
                                        const superposition &prm_superposition, ///< TODOCUMENT
                                        const pdb_list      &prm_pdbs,          ///< TODOCUMENT
                                        const str_vec       &prm_names          ///< TODOCUMENT
                                        ) const {
	const size_t num_pdbs = prm_pdbs.size();
	for (const size_t &pdb_ctr : indices( num_pdbs ) ) {
		prm_os << "cmd.read_pdbstr(\"\"\"";
		ostringstream superposed_pdb_ss;
		write_superposed_pdb_to_ostream( superposed_pdb_ss, prm_superposition, prm_pdbs[pdb_ctr], pdb_ctr );
		prm_os << replace_all_copy(superposed_pdb_ss.str(), "\n", "\\\n");
		prm_os << "\"\"\",\"" << prm_names[pdb_ctr] << "\")\n";
	}
	prm_os << "hide all\n";
	prm_os << "set cartoon_rect_length  = " << pymol_tools::pymol_size( 2, 1.50,  100, 0.090,  num_pdbs ) << "\n";
	prm_os << "set cartoon_rect_width   = " << pymol_tools::pymol_size( 2, 0.40,  100, 0.024,  num_pdbs ) << "\n";
	prm_os << "set cartoon_oval_length  = " << pymol_tools::pymol_size( 2, 1.35,  100, 0.081,  num_pdbs ) << "\n";
	prm_os << "set cartoon_oval_width   = " << pymol_tools::pymol_size( 2, 0.25,  100, 0.015,  num_pdbs ) << "\n";
	prm_os << "set cartoon_loop_radius  = " << pymol_tools::pymol_size( 2, 0.20,  100, 0.036,  num_pdbs ) << "\n";
	prm_os << "set cartoon_helix_radius = " << pymol_tools::pymol_size( 2, 2.00,  100, 0.120,  num_pdbs ) << R"(
set cartoon_cylindrical_helices, 0
set cartoon_fancy_helices, 0
bg_color white
color    black
)";
}

/// \brief TODOCUMENT
void chimera_viewer::do_define_colour(ostream              &prm_os,         ///< TODOCUMENT
                                      const display_colour &prm_colour,    ///< TODOCUMENT
                                      const string         &prm_colour_name ///< TODOCUMENT
                                      ) const {
	prm_os << "set_color "
	       << prm_colour_name
	       << ", ["
	       << comma_separated_string_of_display_colour(prm_colour)
	       << "]\n";
}

/// \brief Get a string for colouring the base (ie everything) in the colour that has previously been defined with the specified name
string chimera_viewer::do_get_colour_base_str(const string &prm_colour_name ///< The previously-defined colour with which to colour the base
                                              ) const {
	return "colour " + prm_colour_name + "\n";
}

/// \brief TODOCUMENT
string chimera_viewer::do_get_colour_pdb_str(const string &prm_colour_name, ///< TODOCUMENT
                                             const string &prm_pdb_name     ///< TODOCUMENT
                                             ) const {
	return "colour "
		+ prm_colour_name
		+ ", "
		+ prm_pdb_name
		+ "\n";
}

/// \brief TODOCUMENT
///
/// This splits the list of residues into batches of RESIDUE_BATCH_SIZE
/// because PyMOL just ignores a list of residues that's too long
string chimera_viewer::do_get_colour_pdb_residues_str(const string         &prm_colour_name, ///< TODOCUMENT
                                                      const string         &prm_pdb_name,    ///< TODOCUMENT
                                                      const residue_id_vec &prm_residue_ids  ///< TODOCUMENT
                                                      ) const {
	const size_t num_res_names   = prm_residue_ids.size();
	const size_t num_res_batches = num_batches( num_res_names, RESIDUE_BATCH_SIZE, broken_batch_tol::PERMIT );
	return "colour " + prm_colour_name + ", "
		+ join(
			transform_build<str_vec>(
				indices( num_res_batches ),
				[&] (const size_t &batch_idx) {
					return pymol_tools::pymol_res_seln_str(
						prm_pdb_name,
						copy_build<residue_id_vec>(
							batch_subrange( prm_residue_ids, RESIDUE_BATCH_SIZE, batch_idx, broken_batch_tol::PERMIT )
						)
					);
				}
			),
			" or "
		)
		+ "\n";
}

/// \brief TODOCUMENT
void chimera_viewer::do_write_alignment_extras(ostream                     &prm_os,                   ///< TODOCUMENT
                                               const superposition_context &prm_superposition_context ///< TODOCUMENT
                                               ) const {
	write_chimera_global_alignment( prm_os, prm_superposition_context );
//	write_chimera_pair_alignments ( prm_os, prm_superposition_context );
}

/// \brief TODOCUMENT
void chimera_viewer::do_write_end(ostream          &prm_os,  ///< TODOCUMENT
                                  const string_ref &prm_name ///< TODOCUMENT
                                  ) const {
	prm_os << R"(show cartoon
set cartoon_smooth_loops,1
show_as sticks, organic
colour black, organic
select protein, ( byres polymer & name CA )
select organic, organic
select nucleic, ( byres polymer & name P  )
deselect
reset
set field_of_view, 25
set label_size, -0.6
set line_width,  6
set dash_width,  0.7
set dash_radius, 0.02
set seq_view_label_mode, 1
set ribbon_width, 1.5
zoom protein

feedback enable,all,output
)"
	<< (
		prm_name.empty()
		? ""
		: ( "print \"" + prm_name.to_string() + "\"\n" )
	);
}

/// \brief TODOCUMENT
string chimera_viewer::parse_residue_id_for_chimera(const residue_id &prm_residue_id ///< TODOCUMENT
                                                    ) {
	return replace_all_copy( to_string( prm_residue_id ), "-", "\\-" );
}

/// \brief TODOCUMENT
str_vec chimera_viewer::parse_residue_ids_for_chimera(const residue_id_vec &prm_residue_ids ///< TODOCUMENT
                                                      ) {
	str_vec new_residue_ids;
	new_residue_ids.reserve( prm_residue_ids.size() );
	for (const residue_id &the_residue_id : prm_residue_ids) {
		new_residue_ids.push_back( parse_residue_id_for_chimera( the_residue_id ) );
	}
	return new_residue_ids;
}

