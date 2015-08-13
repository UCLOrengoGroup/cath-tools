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

#include "chimera_viewer.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "alignment/alignment.h"
#include "common/batch/batch_functions.h"
#include "common/c++14/cbegin_cend.h"
#include "display/display_colour/display_colour.h"
#include "display/viewer/pymol/pymol_tools.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"
#include "superposition/io/superposition_io.h"
#include "superposition/superposition_context.h"

#include <algorithm>
#include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::detail;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sup;
using namespace std;

using boost::adjacency_list;
using boost::algorithm::join;
using boost::algorithm::replace_all_copy;
using boost::edge_weight_t;
using boost::graph_traits;
using boost::kruskal_minimum_spanning_tree;
using boost::lexical_cast;
using boost::no_property;
using boost::property;
using boost::undirectedS;
using boost::vecS;

/// \brief TODOCUMENT
const size_t chimera_viewer::RESIDUE_BATCH_SIZE( 200 );

/// \brief TODOCUMENT
///
/// \relates chimera_viewer
///
/// \todo Refactor out any similarities between write_chimera_pair_alignments() and write_chimera_global_alignment()
///
/// \todo Separate out the code that identifies the residue links in the alignment from the code that writes
///       this information in Chimera-ese (so that the first bit of code can be reused with other viewers)
void cath::detail::write_chimera_pair_alignments(ostream                     &arg_os,                   ///< TODOCUMENT
                                                 const superposition_context &arg_superposition_context ///< TODOCUMENT
                                                 ) {
	if ( ! arg_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot write PyMOL pair alignments for superposition_context with no alignment"));
	}
	const alignment     &the_alignment     = arg_superposition_context.get_alignment_cref();
// 	const superposition &the_superposition = arg_superposition_context.get_superposition_cref();
	const pdb_list      &pdbs              = arg_superposition_context.get_pdbs_cref();
	const str_vec        names             = clean_names_for_viewer( arg_superposition_context );
	
	// Grab some basic details
	const alignment::size_type num_entries   = min( the_alignment.num_entries(), names.size() );
	const alignment::size_type aln_length    = the_alignment.length();
	const residue_name_vec_vec residue_names = get_backbone_complete_residue_names_of_first_chains( pdbs );

	for (size_t entry_ctr_a = 0; entry_ctr_a < num_entries; ++entry_ctr_a) {
		const string           &name_a          = names        [ entry_ctr_a ];
		const residue_name_vec &residue_names_a = residue_names[ entry_ctr_a ];

		for (size_t entry_ctr_b = entry_ctr_a+1; entry_ctr_b < num_entries; ++entry_ctr_b) {
			const string           &name_b          = names        [ entry_ctr_b ];
			const residue_name_vec &residue_names_b = residue_names[ entry_ctr_b ];

			bool added_pair_distances(false);
			for (alignment::size_type aln_posn_ctr = 0; aln_posn_ctr < aln_length; ++aln_posn_ctr) {
				const opt_aln_posn position_a = the_alignment.position_of_entry_of_index( entry_ctr_a, aln_posn_ctr );
				const opt_aln_posn position_b = the_alignment.position_of_entry_of_index( entry_ctr_b, aln_posn_ctr );
				if ( position_a && position_b) {
					added_pair_distances = true;
					if ( *position_a >= residue_names_a.size() ) {
						BOOST_THROW_EXCEPTION(invalid_argument_exception(
							"Whilst adding alignment extras in chimera_viewer, residue index "
							+ lexical_cast<string>( *position_a )
							+ " is out of range "
							+ lexical_cast<string>(residue_names_a.size())
						));
					}
					if ( *position_b >= residue_names_b.size()) {
						BOOST_THROW_EXCEPTION(invalid_argument_exception(
							"Whilst adding alignment extras in chimera_viewer, residue index "
							+ lexical_cast<string>( *position_b )
							+ " is out of range "
							+ lexical_cast<string>(residue_names_b.size())
						));
					}
					const residue_name &residue_name_a = residue_names_a[ *position_a ];
					const residue_name &residue_name_b = residue_names_b[ *position_b ];

					arg_os << "distance ";
					arg_os << name_a;
					arg_os << "_";
					arg_os << name_b;
					arg_os << "_alignment, /";
					arg_os << name_a;
					arg_os << "///";
					arg_os << chimera_viewer::parse_residue_name_for_chimera( residue_name_a );
					arg_os << "/CA/, /";
					arg_os << name_b;
					arg_os << "///";
//					arg_os << residue_name_b;
					arg_os << chimera_viewer::parse_residue_name_for_chimera( residue_name_b );
					arg_os << "/CA/\n";
				}
			}

			if (added_pair_distances) {
				arg_os << "disable ";
				arg_os << name_a;
				arg_os << "_";
				arg_os << name_b;
				arg_os << "_alignment\n";
			}
		}
	}
	arg_os << "hide labels\n";
	arg_os << "set dash_gap,    0.0\n";
	arg_os << "set dash_color,  black\n";
	arg_os << "set dash_radius, 0.05\n";

}

/// \brief TODOCUMENT
///
/// \relates chimera_viewer
///
/// \todo Refactor out any similarities between write_chimera_pair_alignments() and write_chimera_global_alignment()
///
/// \todo Separate out the code that identifies the residue links in each alignment from the code that writes
///       this information in Chimera-ese (so that the first bit of code can be reused with other viewers)
void cath::detail::write_chimera_global_alignment(ostream                     &arg_os,                   ///< TODOCUMENT
                                                  const superposition_context &arg_superposition_context ///< TODOCUMENT
                                                  ) {
	if ( ! arg_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot write PyMOL global alignment for superposition_context with no alignment"));
	}
	const alignment     &the_alignment     = arg_superposition_context.get_alignment_cref();
	const superposition &the_superposition = arg_superposition_context.get_superposition_cref();
	const pdb_list      &pdbs              = arg_superposition_context.get_pdbs_cref();
	const str_vec        names             = clean_names_for_viewer( arg_superposition_context );
	
	// Grab some basic details
	const alignment::size_type num_entries   = min( the_alignment.num_entries(), names.size() );
	const alignment::size_type aln_length    = the_alignment.length();
	const residue_name_vec_vec residue_names = get_backbone_complete_residue_names_of_first_chains( pdbs );

	/// ???
	bool added_distances(false);
	for (alignment::size_type aln_index = 0; aln_index < aln_length; ++aln_index) {
		// Prepare some type aliases that are useful for this
		using Graph     = adjacency_list<vecS, vecS, undirectedS, no_property, property <edge_weight_t, double>>;
		using edge_desc = graph_traits<Graph>::edge_descriptor;

		// Prepare edges and distances to be loaded into the graph
		size_size_pair_vec edges;
		doub_vec           distances;

		// ????
		const size_t num_present_posns = num_present_positions_of_index( the_alignment, aln_index );
		for (size_t entry_a = 0; entry_a < num_entries; ++entry_a) {
			for (size_t entry_b = entry_a + 1; entry_b < num_entries; ++entry_b) {

				const opt_aln_posn posn_a = the_alignment.position_of_entry_of_index( entry_a, aln_index );
				const opt_aln_posn posn_b = the_alignment.position_of_entry_of_index( entry_b, aln_index );

				if ( posn_a && posn_b ) {
					added_distances = true;

//					cerr << "(1) At alignment index\t" << aln_index << ",\tabout to get CA coord from\t" << entry_a << "\t(\t" << names[ entry_a ] << "\t), index:\t" << *posn_a << endl;
//					cerr << "(2) At alignment index\t" << aln_index << ",\tabout to get CA coord from\t" << entry_b << "\t(\t" << names[ entry_b ] << "\t), index:\t" << *posn_b << endl;
					const coord         ca_a        = pdbs[ entry_a ].get_residue_ca_coord_of_backbone_complete_index( *posn_a );
					const coord         ca_b        = pdbs[ entry_b ].get_residue_ca_coord_of_backbone_complete_index( *posn_b );
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
			const residue_name  &res_name_a  = residue_names[ entry_a ][ res_index_a ];
			const residue_name  &res_name_b  = residue_names[ entry_b ][ res_index_b ];

			arg_os << "distance ";
			arg_os << "alignment, /";
			arg_os << name_a;
			arg_os << "///";
			arg_os << chimera_viewer::parse_residue_name_for_chimera( res_name_a );
			arg_os << "/CA/, /";
			arg_os << name_b;
			arg_os << "///";
			arg_os << chimera_viewer::parse_residue_name_for_chimera( res_name_b );
			arg_os << "/CA/\n";
		}
	}
	if (added_distances) {
		arg_os << "disable alignment\n";
	}

	// ????
	if ( the_alignment.is_scored() ) {
		using bool_str_str_vec_map_pair = pair<bool, str_str_vec_map>;
		using bool_str_str_vec_map_map = map <bool, str_str_vec_map>;
		bool_str_str_vec_map_map core_res_names_of_entry_name;
		const alignment_residue_scores &the_scores = the_alignment.get_alignment_residue_scores();
		for (size_t entry = 0; entry < num_entries; ++entry) {
			const string &entry_name = names[ entry ];
			for (alignment::size_type index = 0; index < aln_length; ++index) {
				if ( has_score( the_scores, entry, index ) ) {
					const float_score_type  the_score  = get_score( the_scores, entry, index, true, true );
					const bool              is_core    = ( the_score > 0.25 );
					const aln_posn_type     the_posn   = get_position_of_entry_of_index( the_alignment, entry, index );
					const residue_name     &res_name   = residue_names[ entry ][ the_posn ];
					core_res_names_of_entry_name[ is_core ][ entry_name ].push_back( chimera_viewer::parse_residue_name_for_chimera( res_name ) );
				}
			}
		}
		for (const bool_str_str_vec_map_pair &core_data : core_res_names_of_entry_name) {
			const bool            &is_core                 = core_data.first;
			const string           core_name               = ( is_core ? "core" : "noncore" );
			const str_str_vec_map &res_names_of_entry_name = core_data.second;
			str_vec selection_strings;
			for (const str_str_vec_pair &entry_name_and_res_names : res_names_of_entry_name) {
				const string  &entry_name = entry_name_and_res_names.first;
				const str_vec &res_names  = entry_name_and_res_names.second;

				const size_t num_res_names   = res_names.size();
				const size_t num_res_batches = num_batches( num_res_names, chimera_viewer::RESIDUE_BATCH_SIZE, broken_batch_tol::PERMIT );
				for (size_t batch_ctr = 0; batch_ctr < num_res_batches; ++batch_ctr) {
					string batch_string( "/" + entry_name + "///" );
					const size_size_pair begin_and_end = batch_begin_and_end( num_res_names, chimera_viewer::RESIDUE_BATCH_SIZE, batch_ctr, broken_batch_tol::PERMIT );
					for (size_t res_index = begin_and_end.first; res_index < begin_and_end.second; ++res_index) {
						batch_string += ( res_index > 0 ? "+" : "" );
						batch_string += res_names[ res_index ];
					}
					selection_strings.push_back( batch_string );
				}
//				selection_strings.push_back( "/" + entry_name + "///" + join( res_names,  "+"  ) );

			}
			arg_os << "select " << core_name << ", ( " << join( selection_strings, " or " ) << " )\n";
		}
		arg_os << "deselect\n";
	}

	arg_os << "hide labels\n";
	arg_os << "set dash_gap, 0.0\n";
	arg_os << "set dash_color, black\n";
	arg_os << "set dash_radius, 0.05\n";
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
void chimera_viewer::do_write_start(ostream &arg_os ///< TODOCUMENT
                                    ) const {
	arg_os << "feedback disable,all,output\n";
}

/// \brief TODOCUMENT
void chimera_viewer::do_write_load_pdbs(ostream             &arg_os,            ///< TODOCUMENT
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
	arg_os << "set cartoon_rect_length  = " << pymol_tools::pymol_size( 2, 1.50,  100, 0.090,  num_pdbs ) << "\n";
	arg_os << "set cartoon_rect_width   = " << pymol_tools::pymol_size( 2, 0.40,  100, 0.024,  num_pdbs ) << "\n";
	arg_os << "set cartoon_oval_length  = " << pymol_tools::pymol_size( 2, 1.35,  100, 0.081,  num_pdbs ) << "\n";
	arg_os << "set cartoon_oval_width   = " << pymol_tools::pymol_size( 2, 0.25,  100, 0.015,  num_pdbs ) << "\n";
	arg_os << "set cartoon_loop_radius  = " << pymol_tools::pymol_size( 2, 0.20,  100, 0.036,  num_pdbs ) << "\n";
	arg_os << "set cartoon_helix_radius = " << pymol_tools::pymol_size( 2, 2.00,  100, 0.120,  num_pdbs ) << "\n";
	arg_os << "bg_color white\n";
	arg_os << "color    black\n";
}

/// \brief TODOCUMENT
void chimera_viewer::do_define_colour(ostream              &arg_os,         ///< TODOCUMENT
                                      const display_colour &arg_colour,    ///< TODOCUMENT
                                      const string         &arg_colour_name ///< TODOCUMENT
                                      ) const {
	arg_os << "set_color "
	       << arg_colour_name
	       << ", ["
	       << comma_separated_string_of_display_colour(arg_colour)
	       << "]\n";
}

/// \brief TODOCUMENT
void chimera_viewer::do_colour_pdb(ostream      &arg_os,          ///< TODOCUMENT
                                   const string &arg_colour_name, ///< TODOCUMENT
                                   const string &arg_pdb_name     ///< TODOCUMENT
                                   ) const {
	arg_os << "colour "
	       << arg_colour_name
	       << ", "
	       << arg_pdb_name
	       << "\n";
}

/// \brief TODOCUMENT
///
/// This splits the list of residues into batches of RESIDUE_BATCH_SIZE
/// because PyMOL just ignores a list of residues that's too long
void chimera_viewer::do_colour_pdb_residues(ostream                &arg_os,           ///< TODOCUMENT
                                            const string           &arg_colour_name,  ///< TODOCUMENT
                                            const string           &arg_pdb_name,     ///< TODOCUMENT
                                            const residue_name_vec &arg_residue_names ///< TODOCUMENT
                                            ) const {
	arg_os << "colour " << arg_colour_name << ",";

	const size_t num_res_names   = arg_residue_names.size();
	const size_t num_res_batches = num_batches( num_res_names, RESIDUE_BATCH_SIZE, broken_batch_tol::PERMIT );
	for (size_t batch_ctr = 0; batch_ctr < num_res_batches; ++batch_ctr) {
		arg_os << ( batch_ctr > 0 ? " or " : " " );
		arg_os << "/" << arg_pdb_name << "///";
		const size_size_pair begin_and_end = batch_begin_and_end( num_res_names, RESIDUE_BATCH_SIZE, batch_ctr, broken_batch_tol::PERMIT );
		for (size_t res_index = begin_and_end.first; res_index < begin_and_end.second; ++res_index) {
			arg_os << ( res_index > 0 ? "+" : "" );
			arg_os << arg_residue_names[ res_index ];
		}
		arg_os << "//";
	}
	arg_os << "\n";
}

/// \brief TODOCUMENT
void chimera_viewer::do_write_alignment_extras(ostream                     &arg_os,                   ///< TODOCUMENT
                                               const superposition_context &arg_superposition_context ///< TODOCUMENT
                                               ) const {
	write_chimera_global_alignment( arg_os, arg_superposition_context );
//	write_chimera_pair_alignments ( arg_os, arg_superposition_context );
}

/// \brief TODOCUMENT
void chimera_viewer::do_write_end(ostream &arg_os ///< TODOCUMENT
                                  ) const {
	arg_os << "show cartoon\n";
	arg_os << "set cartoon_smooth_loops,1\n";
	arg_os << "reset\n";
	arg_os << "set field_of_view, 25\n";
	arg_os << "set label_size, -0.6\n";
	arg_os << "set line_width,  6\n";
	arg_os << "set dash_width,  0.7\n";
	arg_os << "set dash_radius, 0.02\n";
	arg_os << "set seq_view_label_mode, 1\n";
	arg_os << "set ribbon_width, 1.5\n";
	arg_os << "orient\n";
	arg_os << "feedback enable,all,output\n";
}

/// \brief TODOCUMENT
string chimera_viewer::parse_residue_name_for_chimera(const residue_name &arg_residue_name ///< TODOCUMENT
                                                      ) {
	return replace_all_copy( lexical_cast<string>( arg_residue_name ), "-", "\\-" );
}

/// \brief TODOCUMENT
str_vec chimera_viewer::parse_residue_names_for_chimera(const residue_name_vec &arg_residue_names ///< TODOCUMENT
                                                        ) {
	str_vec new_residue_names;
	new_residue_names.reserve( arg_residue_names.size() );
	for (const residue_name &the_residue_name : arg_residue_names) {
		new_residue_names.push_back( parse_residue_name_for_chimera( the_residue_name ) );
	}
	return new_residue_names;
}

