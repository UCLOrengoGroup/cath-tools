/// \file
/// \brief The pymol_viewer class definitions

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

#include "pymol_viewer.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "alignment/alignment.hpp"
#include "cath_tools_git_version.hpp"
#include "chopping/region/region.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/batch/batch_functions.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "display/display_colourer/display_colourer.hpp"
#include "display/viewer/pymol/pymol_tools.hpp"
#include "display_colour/display_colour.hpp"
#include "exception/invalid_argument_exception.hpp"
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
using boost::lexical_cast;
using boost::no_property;
using boost::property;
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

constexpr size_t pymol_viewer::RESIDUE_BATCH_SIZE;

/// \brief TODOCUMENT
///
/// \relates pymol_viewer
///
/// \todo Refactor out any similarities between write_pymol_pair_alignments() and write_pymol_global_alignment()
///
/// \todo Separate out the code that identifies the residue links in the alignment from the code that writes
///       this information in Pymol-ese (so that the first bit of code can be reused with other viewers)
void cath::detail::write_pymol_pair_alignments(ostream                     &arg_os,                   ///< TODOCUMENT
                                               const superposition_context &arg_superposition_context ///< TODOCUMENT
                                               ) {
	if ( ! arg_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot write PyMOL pair alignments for superposition_context with no alignment"));
	}
	const alignment     &the_alignment     = arg_superposition_context.get_alignment();
// 	const superposition &the_superposition = arg_superposition_context.get_superposition();
	const pdb_list       pdbs              = get_restricted_pdbs( arg_superposition_context );
	const str_vec        names             = clean_names_for_viewer( arg_superposition_context );
	
	// Grab some basic details
	const alignment::size_type num_entries = min( the_alignment.num_entries(), names.size() );
	const alignment::size_type aln_length  = the_alignment.length();
	const residue_id_vec_vec   residue_ids = get_backbone_complete_residue_ids_of_first_chains( pdbs );

	for (size_t entry_ctr_a = 0; entry_ctr_a < num_entries; ++entry_ctr_a) {
		const string         &name_a        = names      [ entry_ctr_a ];
		const residue_id_vec &residue_ids_a = residue_ids[ entry_ctr_a ];

		for (size_t entry_ctr_b = entry_ctr_a+1; entry_ctr_b < num_entries; ++entry_ctr_b) {
			const string         &name_b        = names      [ entry_ctr_b ];
			const residue_id_vec &residue_ids_b = residue_ids[ entry_ctr_b ];

			bool added_pair_distances(false);
			for (alignment::size_type aln_posn_ctr = 0; aln_posn_ctr < aln_length; ++aln_posn_ctr) {
				const aln_posn_opt position_a = the_alignment.position_of_entry_of_index( entry_ctr_a, aln_posn_ctr );
				const aln_posn_opt position_b = the_alignment.position_of_entry_of_index( entry_ctr_b, aln_posn_ctr );
				if ( position_a && position_b) {
					added_pair_distances = true;
					if ( *position_a >= residue_ids_a.size() ) {
						BOOST_THROW_EXCEPTION(invalid_argument_exception(
							"Whilst adding alignment extras in pymol_viewer, residue index "
							+ lexical_cast<string>( *position_a )
							+ " is out of range "
							+ lexical_cast<string>(residue_ids_a.size())
						));
					}
					if ( *position_b >= residue_ids_b.size()) {
						BOOST_THROW_EXCEPTION(invalid_argument_exception(
							"Whilst adding alignment extras in pymol_viewer, residue index "
							+ lexical_cast<string>( *position_b )
							+ " is out of range "
							+ lexical_cast<string>(residue_ids_b.size())
						));
					}
					const residue_id &residue_id_a = residue_ids_a[ *position_a ];
					const residue_id &residue_id_b = residue_ids_b[ *position_b ];

					arg_os << "distance "
						+ name_a
						+ "_"
						+ name_b
						+ "_alignment, "
						+ pymol_tools::pymol_res_seln_str( name_a, { residue_id_a }, "CA"s )
						+ ", "
						+ pymol_tools::pymol_res_seln_str( name_b, { residue_id_b }, "CA"s )
						+ "\n";
				}
			}

			if (added_pair_distances) {
				arg_os << "disable "
					+ name_a
					+ "_"
					+ name_b
					+ "_alignment\n";
			}
		}
	}
	arg_os << "hide labels\n"
		"set dash_gap,    0.0\n"
		"set dash_color,  black\n"
		"set dash_radius, 0.05\n";

}

/// \brief TODOCUMENT
///
/// \relates pymol_viewer
///
/// \todo Refactor out any similarities between write_pymol_pair_alignments() and write_pymol_global_alignment()
///
/// \todo Separate out the code that identifies the residue links in each alignment from the code that writes
///       this information in Pymol-ese (so that the first bit of code can be reused with other viewers)
void cath::detail::write_pymol_global_alignment(ostream                     &arg_os,                   ///< TODOCUMENT
                                                const superposition_context &arg_superposition_context ///< TODOCUMENT
                                                ) {
	if ( ! arg_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot write PyMOL global alignment for superposition_context with no alignment"));
	}
	const alignment     &the_alignment     = arg_superposition_context.get_alignment();
	const superposition &the_superposition = arg_superposition_context.get_superposition();
	const pdb_list       pdbs              = get_restricted_pdbs( arg_superposition_context );
	const str_vec        names             = clean_names_for_viewer( arg_superposition_context );
	
	// Grab some basic details
	const alignment::size_type num_entries = min( the_alignment.num_entries(), names.size() );
	const alignment::size_type aln_length  = the_alignment.length();
	const residue_id_vec_vec   residue_ids = get_backbone_complete_residue_ids( pdbs );

	/// ???
	bool added_distances(false);
	for (alignment::size_type aln_index = 0; aln_index < aln_length; ++aln_index) {
		// Prepare some type aliases that are useful for this
		using Graph = adjacency_list < vecS, vecS, undirectedS, no_property, property <edge_weight_t, double> >;
		using edge_desc = graph_traits < Graph >::edge_descriptor;

		// Prepare edges and distances to be loaded into the graph
		size_size_pair_vec edges;
		doub_vec           distances;

		// ????
		const size_t num_present_posns = num_present_positions_of_index( the_alignment, aln_index );
		for (size_t entry_a = 0; entry_a < num_entries; ++entry_a) {
			for (size_t entry_b = entry_a + 1; entry_b < num_entries; ++entry_b) {

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

			arg_os << "distance alignment, "
				+ pymol_tools::pymol_res_seln_str( name_a, { res_name_a }, "CA"s )
				+ ", "
				+ pymol_tools::pymol_res_seln_str( name_b, { res_name_b }, "CA"s )
				+ "\n";
		}
	}
	if (added_distances) {
		arg_os << "disable alignment\n";
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
		for (size_t entry = 0; entry < num_entries; ++entry) {
			const string &entry_name = names[ entry ];
			for (alignment::size_type index = 0; index < aln_length; ++index) {
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
		arg_os << join(
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
						const size_t num_res_batches = num_batches( num_res_ids, pymol_viewer::RESIDUE_BATCH_SIZE, broken_batch_tol::PERMIT );
						for (const size_t &batch_ctr : indices( num_res_batches ) ) {
							const size_size_pair begin_and_end = batch_begin_and_end(
								num_res_ids,
								pymol_viewer::RESIDUE_BATCH_SIZE,
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
		arg_os << "deselect\n";
	}

	arg_os << "hide labels\n";
	arg_os << "set dash_gap, 0.0\n";
	arg_os << "set dash_color, black\n";
	arg_os << "set dash_radius, 0.05\n";
}

/// \brief Write PyMOL commands to store the current scene under the specified label to the specified ostream
void pymol_viewer::record_scene(ostream           &arg_os,             ///< The ostream to which to write the PyMOL commands
                                const std::string &arg_colouring_label ///< The label under which the scene should be stored
                                ) {
	arg_os
		<< "scene F"
		<< scene_count
		<< ", store, message=\""
		<< arg_colouring_label
		<< "\", color=1, view=0, active=0, rep=0, frame=0\n";
	++scene_count;
}

/// \brief TODOCUMENT
string pymol_viewer::do_default_executable() const {
	return "pymol";
}

/// \brief TODOCUMENT
string pymol_viewer::do_default_file_extension() const {
	return ".pml";
}

/// \brief TODOCUMENT
void pymol_viewer::do_write_start(ostream &arg_os ///< TODOCUMENT
                                  ) const {
	arg_os << "feedback disable,all,output\n";
}

/// \brief TODOCUMENT
void pymol_viewer::do_write_load_pdbs(ostream             &arg_os,            ///< TODOCUMENT
                                      const superposition &arg_superposition, ///< TODOCUMENT
                                      const pdb_list      &arg_pdbs,          ///< TODOCUMENT
                                      const str_vec       &arg_names          ///< TODOCUMENT
                                      ) const {
	const size_t num_pdbs = arg_pdbs.size();
	for (size_t pdb_ctr = 0; pdb_ctr < num_pdbs; ++pdb_ctr) {
		arg_os << "cmd.read_pdbstr(\"\"\"";
		ostringstream superposed_pdb_ss;
		write_superposed_pdb_to_ostream( superposed_pdb_ss, arg_superposition, arg_pdbs[pdb_ctr], pdb_ctr );
		arg_os << replace_all_copy(superposed_pdb_ss.str(), "\n", "\\\n");
		arg_os << "\"\"\",\"" << arg_names[pdb_ctr] << "\")\n";
	}
	arg_os << "hide all\n";
	arg_os << "set cartoon_rect_length  = " << pymol_tools::pymol_size(2, 1.50,  100, 0.090,  num_pdbs) << "\n";
	arg_os << "set cartoon_rect_width   = " << pymol_tools::pymol_size(2, 0.40,  100, 0.024,  num_pdbs) << "\n";
	arg_os << "set cartoon_oval_length  = " << pymol_tools::pymol_size(2, 1.35,  100, 0.081,  num_pdbs) << "\n";
	arg_os << "set cartoon_oval_width   = " << pymol_tools::pymol_size(2, 0.25,  100, 0.015,  num_pdbs) << "\n";
	arg_os << "set cartoon_loop_radius  = " << pymol_tools::pymol_size(2, 0.20,  100, 0.036,  num_pdbs) << "\n";
	arg_os << "set cartoon_helix_radius = " << pymol_tools::pymol_size(2, 2.00,  100, 0.120,  num_pdbs) << "\n";
	arg_os << "bg_color white\n";
	arg_os << "color    black\n";
}

/// \brief TODOCUMENT
void pymol_viewer::do_define_colour(ostream              &arg_os,         ///< TODOCUMENT
                                    const display_colour &arg_colour,     ///< TODOCUMENT
                                    const string         &arg_colour_name ///< TODOCUMENT
                                    ) const {
	arg_os << "set_color "
		+ arg_colour_name
		+ ", ["
		+ comma_separated_string_of_display_colour( arg_colour )
		+ "]\n";
}

/// \brief Specify that pymol_view does accept multiple colourings (because it can store them as scenes)
bool pymol_viewer::do_accepts_multiple_colourings() const {
	return true;
}

/// \brief Write PyMOL commands to the specified ostream to prepare for a new, specified colouring
void pymol_viewer::do_begin_colouring(ostream                &arg_os,          ///< The ostream to which the PyMOL commands should be written
                                      const display_colourer &/*arg_colourer*/ ///< The display_colourer to be used for the colouring that is beginning
                                      ) {
	// If this is the first colouring, precede it by colouring by secondary structure and storing that as the first scene
	if ( scene_count == 1 ) {
		arg_os << R"(color black
color density, ss s
color rutherfordium, ss h
)";
		record_scene( arg_os, "Colour by secondary structure" );
	}
}

/// \brief Get a string for colouring the base (ie everything) in the colour that has previously been defined with the specified name
string pymol_viewer::do_get_colour_base_str(const string &arg_colour_name ///< The previously-defined colour with which to colour the base
                                            ) const {
	return "colour " + arg_colour_name + "\n";
}

/// \brief TODOCUMENT
string pymol_viewer::do_get_colour_pdb_str(const string &arg_colour_name, ///< TODOCUMENT
                                           const string &arg_pdb_name     ///< TODOCUMENT
                                           ) const {
	return "colour "
		+ arg_colour_name
		+ R"(, ")"
		+ arg_pdb_name
		+ "\"\n";
}

/// \brief TODOCUMENT
///
/// This splits the list of residues into batches of RESIDUE_BATCH_SIZE
/// because PyMOL just ignores a list of residues that's too long
string pymol_viewer::do_get_colour_pdb_residues_str(const string         &arg_colour_name, ///< TODOCUMENT
                                                    const string         &arg_pdb_name,    ///< TODOCUMENT
                                                    const residue_id_vec &arg_residue_ids  ///< TODOCUMENT
                                                    ) const {
	const size_t num_res_names   = arg_residue_ids.size();
	const size_t num_res_batches = num_batches( num_res_names, RESIDUE_BATCH_SIZE, broken_batch_tol::PERMIT );
	return "colour " + arg_colour_name + ", "
		+ join(
			transform_build<str_vec>(
				indices( num_res_batches ),
				[&] (const size_t &batch_idx) {
					return pymol_tools::pymol_res_seln_str(
						arg_pdb_name,
						copy_build<residue_id_vec>(
							batch_subrange( arg_residue_ids, RESIDUE_BATCH_SIZE, batch_idx, broken_batch_tol::PERMIT )
						)
					);
				}
			),
			" or "
		)
		+ "\n";
}

/// \brief Write PyMOL commands to the specified ostream to finish the specified colouring
void pymol_viewer::do_end_colouring(ostream                &arg_os,      ///< The ostream to which the PyMOL commands should be written
                                    const display_colourer &arg_colourer ///< The display_colourer to be used for the colouring that is ending
                                    ) {
	// Store the scene
	record_scene( arg_os, arg_colourer.get_label() );
}

/// \brief TODOCUMENT
void pymol_viewer::do_write_alignment_extras(ostream                     &arg_os,                   ///< TODOCUMENT
                                             const superposition_context &arg_superposition_context ///< TODOCUMENT
                                             ) const {
	write_pymol_global_alignment( arg_os, arg_superposition_context );
//	write_pymol_pair_alignments ( arg_os, arg_superposition_context );
}

/// \brief TODOCUMENT
void pymol_viewer::do_write_end(ostream &arg_os ///< TODOCUMENT
                                ) const {
	const string CATH_TOOLS_VERSION{ CATH_TOOLS_GIT_VERSION };
	const string advert_msg = "Superposition generated by cath-superpose (" + CATH_TOOLS_VERSION + "), one of the cath-tools."
		+ (
			( scene_count > 2 )
				? ( " Use functions keys F1 - F" + std::to_string( scene_count - 1 ) + " to switch between colouring schemes." )
				: ""
		);
	arg_os << R"(scene F)" << scene_count << R"(, store, message="Colour me badd", color=0, view=0, active=0, rep=0, frame=0
show cartoon
set cartoon_smooth_loops,1
show_as sticks, organic
colour black, organic
select organic, organic
deselect
reset
set field_of_view, 25
set label_size, -0.6
set line_width,  6
set dash_width,  0.7
set dash_radius, 0.02
set seq_view_label_mode, 1
set ribbon_width, 1.5
orient ( polymer and not resn A+C+G+T+U )
cmd.wizard( "message", ")" << advert_msg << R"(" );
feedback enable,all,output
)";
}

