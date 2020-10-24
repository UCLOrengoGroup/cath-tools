/// \file
/// \brief The aligned_pair_score_list_factory class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_ALIGNED_PAIR_SCORE_LIST_FACTORY_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_ALIGNED_PAIR_SCORE_LIST_FACTORY_HPP

#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/mpl/vector.hpp>

namespace cath {
	namespace score {
		class aligned_pair_score;
		class aligned_pair_score_list;

		class drmsd_score;
		class gsas_score;
		class length_score;
		class lddt_score;
		class mi_score;
		class rmsd_score;
		class sas_score;
		class sequence_similarity_score;
		class si_score;

		using all_aligned_pair_score_types = boost::mpl::vector< drmsd_score,
		                                                         gsas_score,
		                                                         length_score,
		                                                         lddt_score,
		                                                         mi_score,
		                                                         rmsd_score,
		                                                         sas_score,
		                                                         sequence_similarity_score,
		                                                         si_score >;

		aligned_pair_score_list make_full_aligned_pair_score_list();
		aligned_pair_score_list make_seq_sim_flavours_aligned_pair_score_list();
		aligned_pair_score_list make_ssap_flavours_aligned_pair_score_list();
		aligned_pair_score_list make_one_of_each_type_aligned_pair_score_list();
		aligned_pair_score_list make_default_aligned_pair_score_list();
		aligned_pair_score_list make_old_ssap_aligned_pair_score_list();
	} // namespace score

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_ALIGNED_PAIR_SCORE_LIST_FACTORY_HPP
