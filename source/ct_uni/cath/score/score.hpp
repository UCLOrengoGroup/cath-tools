/// \file
/// \brief Documentation for the cath::score namespace

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_HPP

/// \namespace cath::score
///
/// \brief Tools for measures of similarity between two aligned structures
///
/// Contents
/// --------
///  * \ref cath_score_aligned_pair_score_group, containing:
///    * gsas_score
///    * lddt_score
///  * \ref cath_score_aligned_pair_score_list, containing:
///

/// \defgroup cath_score_aligned_pair_score_group Aligned Pair Scores
/// \brief Specific scores (eg rmsd_score) that calculate some measure of similarity between two proteins and an alignment

// \defgroup cath_score_aligned_pair_score_list Aligned Pair Score List
// \brief Specific scores (eg rmsd_score) that calculate some measure of similarity between two proteins and an alignment

// @{
// aligned_pair_score
// drmsd_score
// gsas_score
// lddt_score
// length_score
// mi_score
// rmsd_score
// sas_score
// si_score
// @}
//
// @{
// aligned_pair_score_list_factory
// aligned_pair_score_list
// aligned_pair_score_value_list
// @}
//
// @{
// score_value_list_outputter/score_value_list_json_outputter
// @}
//
// @{
// length_getter
// length_of_first_getter
// length_of_longer_getter
// length_of_second_getter
// length_of_shorter_getter
// num_aligned_length_getter
// protein_only_length_getter
// sym_protein_only_length_getter
// @}
//
// \subsection cath_scores_implementations Implementations
//
// score/aligned_pair_score/detail/score_common_coord_handler.hpp
// score/aligned_pair_score_list/detail/aligned_pair_score_list_append.hpp
// score/aligned_pair_score_list/score_value_list_outputter/score_value_list_json_outputter.hpp

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_HPP
