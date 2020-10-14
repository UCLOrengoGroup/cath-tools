/// \file
/// \brief The populate_matrix_scan_action class definitions

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

#include "populate_matrix_scan_action.hpp"

#include <boost/numeric/conversion/cast.hpp>

#include "alignment/dyn_prog_align/detail/matrix_plotter/gnuplot_matrix_plotter.hpp"
#include "alignment/dyn_prog_align/dyn_prog_score_source/new_matrix_dyn_prog_score_source.hpp"
#include "common/exception/invalid_argument_exception.hpp"

using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::common;
using namespace cath::scan;
// using namespace std;

using boost::filesystem::path;
using boost::numeric_cast;

/// \brief TODOCUMENT
populate_matrix_scan_action::populate_matrix_scan_action(const index_type &prm_num_residues_a, ///< TODOCUMENT
                                                         const index_type &prm_num_residues_b, ///< TODOCUMENT
                                                         const index_type &prm_structure_a,    ///< TODOCUMENT
                                                         const index_type &prm_structure_b     ///< TODOCUMENT
                                                         ) : the_matrix  ( prm_num_residues_a, doub_vec( prm_num_residues_b, 0.0 ) ),
                                                             structure_a ( prm_structure_a ),
                                                             structure_b ( prm_structure_b ) {
	if ( prm_num_residues_a == 0 && prm_num_residues_b == 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to create populate_matrix_scan_action with 0 residues in either dimension"));
	}
}

//new_matrix_dyn_prog_score_source

/// \brief TODOCUMENT
void populate_matrix_scan_action::plot_to_file(const path     &prm_out_filename, ///< TODOCUMENT
                                               matrix_plotter &prm_plotter       ///< TODOCUMENT
                                               ) const {
	new_matrix_dyn_prog_score_source the_source{ the_matrix, get_length_a(), get_length_b() };
	prm_plotter.plot_scores( the_source );
	prm_plotter.finish( prm_out_filename );
//	return ;
}

/// \brief TODOCUMENT
index_type populate_matrix_scan_action::get_length_a() const {
	return numeric_cast<index_type>( the_matrix.size() );
}

/// \brief TODOCUMENT
index_type populate_matrix_scan_action::get_length_b() const {
	return numeric_cast<index_type>( the_matrix.front().size() );
}

//			new_matrix_dyn_prog_score_source

/// \brief TODOCUMENT
void cath::scan::gnuplot_to_file(const populate_matrix_scan_action &prm_action,      ///< TODOCUMENT
                                 const path                        &prm_out_filename ///< TODOCUMENT
                                 ) {
	gnuplot_matrix_plotter the_plotter{
		prm_action.get_length_a(),
		prm_action.get_length_b()
	};
	prm_action.plot_to_file(
		prm_out_filename,
		the_plotter
	);
}
