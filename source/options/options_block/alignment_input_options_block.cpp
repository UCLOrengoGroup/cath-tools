/// \file
/// \brief The alignment_input_options_block class definitions

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

#include "alignment_input_options_block.h"

#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/optional.hpp>

#include "common/clone/make_uptr_clone.h"
#include "options/acquirer/alignment_acquirer/cora_aln_file_alignment_acquirer.h"
#include "options/acquirer/alignment_acquirer/fasta_aln_file_alignment_acquirer.h"
#include "options/acquirer/alignment_acquirer/residue_name_alignment_acquirer.h"
#include "options/acquirer/alignment_acquirer/ssap_aln_file_alignment_acquirer.h"
#include "options/acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer.h"

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

using boost::assign::ptr_push_back;
using boost::none;
using boost::ptr_vector;

const string alignment_input_options_block::PO_RES_NAME_ALIGN    ( "res-name-align"     );
const string alignment_input_options_block::PO_FASTA_ALIGN_INFILE( "fasta-aln-infile"    );
const string alignment_input_options_block::PO_SSAP_ALIGN_INFILE ( "ssap-aln-infile"    );
const string alignment_input_options_block::PO_CORA_ALIGN_INFILE ( "cora-aln-infile"    );
const string alignment_input_options_block::PO_SSAP_SCORE_INFILE ( "ssap-scores-infile" );

/// \brief A standard do_clone method.
unique_ptr<options_block> alignment_input_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string alignment_input_options_block::do_get_block_name() const {
	return "Alignment source";
}

/// \brief Add this block's options to the provided options_description
void alignment_input_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                  ) {
	arg_desc.add_options()
		( PO_RES_NAME_ALIGN.c_str(),     bool_switch( &residue_name_align   )->default_value(false), "Align residues by simply matching their names (numbers+insert)\n(for multiple models of the same structure)" )
		( PO_FASTA_ALIGN_INFILE.c_str(), value<path>( &fasta_alignment_file ),                       "Read FASTA alignment from file arg"                                                  )
		( PO_SSAP_ALIGN_INFILE.c_str(),  value<path>( &ssap_alignment_file  ),                       "Read SSAP alignment from file arg"                                                   )
		( PO_CORA_ALIGN_INFILE.c_str(),  value<path>( &cora_alignment_file  ),                       "Read CORA alignment from file arg"                                                   )
		( PO_SSAP_SCORE_INFILE.c_str(),  value<path>( &ssap_scores_file     ),                       "Read SSAP scores from file arg\nAssumes all .list alignment files in same directory" );
}

/// \brief TODOCUMENT
opt_str alignment_input_options_block::do_invalid_string() const {
	if ( ! get_fasta_alignment_file().empty() && ! is_acceptable_input_file( get_fasta_alignment_file()    ) ) {
		return "FASTA alignment file " + get_ssap_alignment_file().string() + " is not a valid input file";
	}
	if ( ! get_ssap_alignment_file().empty()  && ! is_acceptable_input_file( get_ssap_alignment_file()    ) ) {
		return "SSAP alignment file " + get_ssap_alignment_file().string() + " is not a valid input file";
	}
	if ( ! get_cora_alignment_file().empty()  && ! is_acceptable_input_file( get_cora_alignment_file()    ) ) {
		return "CORA alignment file " + get_cora_alignment_file().string() + " is not a valid input file";
	}
	if ( ! get_ssap_scores_file().empty()     && ! is_acceptable_input_file( get_ssap_scores_file(), true ) ) {
		return "SSAP scores file "    + get_ssap_scores_file().string()    + " is not a valid input file";
	}
	return none;
}

/// \brief TODOCUMENT
path alignment_input_options_block::get_fasta_alignment_file() const {
	return fasta_alignment_file;
}

/// \brief TODOCUMENT
path alignment_input_options_block::get_ssap_alignment_file() const {
	return ssap_alignment_file;
}

/// \brief TODOCUMENT
path alignment_input_options_block::get_cora_alignment_file() const {
	return cora_alignment_file;
}

/// \brief TODOCUMENT
///
/// Temporarily public so that cath_superpose_options can flick to common_residue_select_all_policy
/// when performing residue name aligning
bool alignment_input_options_block::get_residue_name_align() const {
	return residue_name_align;
}

/// \brief TODOCUMENT
path alignment_input_options_block::get_ssap_scores_file() const {
	return ssap_scores_file;
}

/// \brief Construct suitable alignment_acquirer objects
ptr_vector<alignment_acquirer> alignment_input_options_block::get_alignment_acquirers() const {
	ptr_vector<alignment_acquirer> alignment_acquirers;

	//      If the alignment is to be created by reading a legacy CORA alignment file, do that
	// Else if the alignment is to be created by aligning using residue names,         do that
	// Else if the alignment is to be created by reading a FASTA alignment file,       do that
	// Else if the alignment is to be created by reading a legacy SSAP alignment file, do that
	// Else if the alignment is to be created by reading a list of SSAP scores,        do that
	if ( ! get_cora_alignment_file().empty()  ) {
		ptr_push_back< cora_aln_file_alignment_acquirer    >( alignment_acquirers )( get_cora_alignment_file()  );
	}
	if (   get_residue_name_align()           ) {
		ptr_push_back< residue_name_alignment_acquirer     >( alignment_acquirers )(                            );
	}
	if ( ! get_fasta_alignment_file().empty()  ) {
		ptr_push_back< fasta_aln_file_alignment_acquirer   >( alignment_acquirers )( get_fasta_alignment_file() );
	}
	if ( ! get_ssap_alignment_file().empty()  ) {
		ptr_push_back< ssap_aln_file_alignment_acquirer    >( alignment_acquirers )( get_ssap_alignment_file()  );
	}
	if ( ! get_ssap_scores_file().empty()     ) {
		ptr_push_back< ssap_scores_file_alignment_acquirer >( alignment_acquirers )( get_ssap_scores_file()     );
	}

	return alignment_acquirers;
}
