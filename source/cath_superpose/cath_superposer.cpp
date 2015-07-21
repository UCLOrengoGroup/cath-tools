/// \file
/// \brief The cath_superposer class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "cath_superposer.h"

#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/alignment_context.h"
#include "options/executable/cath_superpose_options/cath_superpose_options.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "options/outputter/alignment_outputter/alignment_outputter.h"
#include "options/outputter/alignment_outputter/alignment_outputter_list.h"
#include "options/outputter/superposition_outputter/superposition_outputter.h"
#include "options/outputter/superposition_outputter/superposition_outputter_list.h"
#include "superposition/superposition_context.h"

using namespace cath;
using namespace cath::align;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

using boost::ptr_vector;

/// \brief Perform a cath-superpose job as specified by the cath_superpose_options argument
///
/// The input and output stream parameters default to cin and cout respectively but are configurable,
/// primarily for testing purposes
void cath_superposer::superpose(const cath_superpose_options &arg_cath_superpose_options, ///< The details of the cath_superpose job to perform
                                istream                      &arg_istream,                ///< The istream from which any stdin-like input should be read
                                ostream                      &arg_stdout,                 ///< The ostream to which any stdout-like output should be written
                                ostream                      &arg_stderr                  ///< The ostream to which any stderr-like output should be written
                                ) {
	// If the options are invalid or specify to do_nothing, then just return
	const string error_or_help_string = arg_cath_superpose_options.get_error_or_help_string();
	if (!error_or_help_string.empty()) {
		arg_stderr << error_or_help_string << endl;
		return;
	}

	// Get the superposition (and context) as specified by the cath_superpose_options
	const superposition_context the_sup_context = arg_cath_superpose_options.get_superposition_context( arg_istream, arg_stderr );

	// If there is an alignment, then for each of the alignment_outputters specified by the cath_superpose_options, output it
	if ( the_sup_context.has_alignment() ) {
		const pdb_list                 &pdbs           = the_sup_context.get_pdbs_cref();
		const str_vec                  &names          = the_sup_context.get_names_cref();
		const alignment                &the_aln        = the_sup_context.get_alignment_cref();
		const alignment_outputter_list  aln_outputters = arg_cath_superpose_options.get_alignment_outputters();
		use_all_alignment_outputters( aln_outputters, alignment_context( pdbs, names, the_aln ), arg_stdout, arg_stderr );
	}

	// For each of the superposition_outputters specified by the cath_superpose_options, output the superposition
	const superposition_outputter_list outputters = arg_cath_superpose_options.get_superposition_outputters();
	use_all_superposition_outputters( outputters, the_sup_context, arg_stdout, arg_stderr );
}
