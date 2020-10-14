/// \file
/// \brief The sas_score class definitions

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

#include "sas_score.hpp"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/serialization/export.hpp>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/less_than_helper.hpp"

#include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::tribool;

BOOST_CLASS_EXPORT(sas_score)

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> sas_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that a lower SAS is generally better
tribool sas_score::do_higher_is_better() const {
	return false;
}

/// \brief Concrete implementation for calculating the SAS of an alignment
score_value sas_score::do_calculate(const alignment &prm_alignment, ///< The pair alignment to be scored
                                    const protein   &prm_protein_a, ///< The protein associated with the first  half of the alignment
                                    const protein   &prm_protein_b  ///< The protein associated with the second half of the alignment
                                    ) const {
	return 100.0 * rmsd.calculate( prm_alignment, prm_protein_a, prm_protein_b )
	             / num_aligned_residues.calculate( prm_alignment, prm_protein_a, prm_protein_b );
}

/// \brief Concrete implementation that describes what this score means
string sas_score::do_description() const {
	return "The RMSD"
	       + rmsd.description_brackets_string()
	       + " times by 100 and divided by the number of aligned residues"
	       + num_aligned_residues.description_brackets_string();
}

/// \brief TODOCUMENT
string sas_score::do_id_name() const {
	return "SAS";
}

/// \brief TODOCUMENT
str_bool_pair_vec sas_score::do_short_name_suffixes() const {
	return rmsd.short_name_suffixes();
}

/// \brief Concrete implementation providing long name
string sas_score::do_long_name() const {
	return "Structural Alignment Score"
	       + rmsd.long_suffix_string()
	       + ". This is measured in angstroms.";
}

/// \brief Concrete implementation for providing a reference to a publication describing this score
string sas_score::do_reference() const {
	return { "Subbiah S, Laurents DV, Levitt M. "
	         "Structural similarity of DNA-binding domains of bacteriophage repressors and the globin core. "
	         "Curr Biol. 1993 Mar;3(3):141-8. PubMed PMID: 15335781." };
}

///// \brief Build an aligned_pair_score of this concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> sas_score::do_build_from_short_name_spec(const string &prm_short_name_spec ///< The short_name_spec that defines any properties that the resulting aligned_pair_score should have
//                                                                        ) const {
//	cerr << "Should build a sas_score from string \"" << prm_short_name_spec << "\"" << endl;
//	return clone();
//}

/// \brief TODOCUMENT
bool sas_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &prm_aligned_pair_score ///< TODOCUMENT
                                                    ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( prm_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief Ctor for sas_score that allows the caller to specify the common_residue_selection_policy and common_atom_selection_policy
sas_score::sas_score(const common_residue_selection_policy &prm_comm_res_seln_pol, ///< The policy to use for selecting common residues
                     const common_atom_selection_policy    &prm_comm_atom_seln_pol  ///< The policy to use for selecting common atoms
                     ) : rmsd                ( prm_comm_res_seln_pol, prm_comm_atom_seln_pol ),
                         num_aligned_residues( num_aligned_length_getter( prm_comm_res_seln_pol ) ) {
}

/// \brief TODOCUMENT
const rmsd_score & sas_score::get_rmsd_score() const {
	return rmsd;
}

/// \brief TODOCUMENT
const length_score & sas_score::get_num_aligned_residues() const {
	return num_aligned_residues;
}

/// \brief TODOCUMENT
///
/// \relates sas_score
bool cath::score::operator<(const sas_score &prm_sas_score_a, ///< TODOCUMENT
                            const sas_score &prm_sas_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( prm_sas_score_a, prm_sas_score_b );
	the_helper.register_comparison_field( &sas_score::get_rmsd_score           );
	the_helper.register_comparison_field( &sas_score::get_num_aligned_residues );
	return final_less_than_result( the_helper );
}

