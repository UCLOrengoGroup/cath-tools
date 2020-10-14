/// \file
/// \brief The gsas_score class definitions

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

#include "gsas_score.hpp"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/serialization/export.hpp>

#include "alignment/alignment.hpp"
#include "alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "alignment/gap/alignment_gap.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/less_than_helper.hpp"

// #include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::align;
using namespace cath::align::gap;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::numeric_cast;
using boost::tribool;

BOOST_CLASS_EXPORT(gsas_score)

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> gsas_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that a lower SAS is generally better
tribool gsas_score::do_higher_is_better() const {
	return false;
}

/// \brief Concrete implementation for calculating the SAS of an alignment
score_value gsas_score::do_calculate(const alignment &prm_alignment, ///< The pair alignment to be scored
                                     const protein   &prm_protein_a, ///< The protein associated with the first  half of the alignment
                                     const protein   &prm_protein_b  ///< The protein associated with the second half of the alignment
                                     ) const {
	const score_value bad_value   = 99.9;

	// Grab the number of gaps and the number of aligned residues
//	const score_value num_gaps    = gap_count_of_alignment( prm_alignment );
	const score_value num_gaps    = numeric_cast<score_value>( get_naive_num_gaps( prm_alignment ) );
	const score_value num_aligned = num_aligned_residues.calculate( prm_alignment, prm_protein_a, prm_protein_b );

	// If the number of gaps meets or exceeds the number of aligned residues, return a bad value (99.9)
	if ( num_gaps >= num_aligned ) {
		return bad_value;
	}

	// Otherwise, return the standard GSAS formula: 100 * rmsd / (num_aligned - num_gaps)
	const score_value rmsd_val = rmsd.calculate( prm_alignment, prm_protein_a, prm_protein_b );
	return 100.0 * rmsd_val / ( num_aligned - num_gaps );
}

/// \brief Concrete implementation that describes what this score means
string gsas_score::do_description() const {
	return "When the number of gaps matches or exceeds the number of aligned residues, 99.9, otherwise, the RMSD"
	       + rmsd.description_brackets_string()
	       + " times by 100 and divided by (the number of aligned residues"
	       + num_aligned_residues.description_brackets_string()
	       + " minus the total number of gaps in the alignment)";
}

/// \brief TODOCUMENT
string gsas_score::do_id_name() const {
	return "GSAS";
}

/// \brief TODOCUMENT
str_bool_pair_vec gsas_score::do_short_name_suffixes() const {
	return rmsd.short_name_suffixes();
}

/// \brief Concrete implementation providing long name
string gsas_score::do_long_name() const {
	return "Gapped Structural Alignment Score"
	       + rmsd.long_suffix_string()
	       + ". This is measured in angstroms.";
}

/// \brief Concrete implementation for providing a reference to a publication describing this score
string gsas_score::do_reference() const {
	return { "Kolodny R, Koehl P, Levitt M. "
	         "Comprehensive evaluation of protein structure alignment methods: scoring by geometric measures. "
	         "J Mol Biol. 2005 Mar 4;346(4):1173-88. Epub 2005 Jan 16. PubMed PMID: 15701525; PubMed Central PMCID: PMC2692023." };
}

///// \brief Build an aligned_pair_score of this concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> gsas_score::do_build_from_short_name_spec(const string &prm_short_name_spec ///< The short_name_spec that defines any properties that the resulting aligned_pair_score should have
//                                                                         ) const {
//	cerr << "Should build a gsas_score from string \"" << prm_short_name_spec << "\"" << endl;
//	return clone();
//}

/// \brief TODOCUMENT
bool gsas_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &prm_aligned_pair_score ///< TODOCUMENT
                                                     ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( prm_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief Ctor for gsas_score that allows the caller to specify the common_residue_selection_policy and common_atom_selection_policy
gsas_score::gsas_score(const common_residue_selection_policy &prm_comm_res_seln_pol, ///< The policy to use for selecting common residues
                       const common_atom_selection_policy    &prm_comm_atom_seln_pol  ///< The policy to use for selecting common atoms
                       ) : rmsd                ( prm_comm_res_seln_pol, prm_comm_atom_seln_pol      ),
                           num_aligned_residues( num_aligned_length_getter( prm_comm_res_seln_pol ) ) {
}

/// \brief TODOCUMENT
const rmsd_score & gsas_score::get_rmsd_score() const {
	return rmsd;
}

/// \brief TODOCUMENT
const length_score & gsas_score::get_num_aligned_residues() const {
	return num_aligned_residues;
}

/// \brief TODOCUMENT
///
/// \relates gsas_score
bool cath::score::operator<(const gsas_score &prm_gsas_score_a, ///< TODOCUMENT
                            const gsas_score &prm_gsas_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( prm_gsas_score_a, prm_gsas_score_b );
	the_helper.register_comparison_field( &gsas_score::get_rmsd_score           );
	the_helper.register_comparison_field( &gsas_score::get_num_aligned_residues );
	return final_less_than_result( the_helper );
}
