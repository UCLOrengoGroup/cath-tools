/// \file
/// \brief The hmmer_hmmsearch_domtblout class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_HMMER_HMMSEARCH_DOMTBLOUT_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_HMMER_HMMSEARCH_DOMTBLOUT_H

#include <boost/filesystem/path.hpp>

#include "common/type_aliases.h"
#include "resolve_hits/resolve_hits_type_aliases.h"
#include "resolve_hits/file/cath_id_score_category.h"

namespace cath { namespace rslv { class read_and_process_mgr; } }

namespace cath {
	namespace rslv {

		// /// \brief Calculate the Jon's scoring scheme for the specified scores from HMMer (and optionally cath_id_score_category)
		// ///
		// /// Issues:
		// ///  * Dividing by a constant should affect the algorithm. So can it be dropped?
		// ///  * Some of the scores in Jon's example file are negative eg, the following two lines are the worst with score -4.0:
		// ///     sp|Q14807|KIF22_HUMAN -            665 3wu2U01_round_3      -             56     2e-15   44.2   0.0   1   2         1         2   -4.0   0.0    13    27   429   443   427   445 0.66 Kinesin-like protein KIF22 OS=Homo sapiens GN=KIF22 PE=1 SV=5
		// ///     sp|Q14807|KIF22_HUMAN -            665 2h58A00_round_1      -            284   4.2e-96  308.5   0.0   2   2      0.36      0.72   -4.0   0.0    83   111   577   605   569   623 0.65 Kinesin-like protein KIF22 OS=Homo sapiens GN=KIF22 PE=1 SV=5
		// ///                      0.00048828125
		// ///    For now, will add 2^-10 = 1/1024 ~= 0.0009765625 which should make scores positive
		// inline constexpr resscr_t jon_score_of_hmmer_scores(const resscr_t               &arg_bitscore,                                      ///< The HMMer bit-score
		//                                                     const resscr_t               &arg_cond_evalue,                                   ///< The HMMer conditional evalue
		//                                                     const resscr_t               &arg_indp_evalue,                                   ///< The HMMer independent evalue
		//                                                     const cath_id_score_category &arg_cath_category = cath_id_score_category::NORMAL ///< (optional) The type of id, allowing for Jon's special handling of certain match ID types
		//                                                     ) {
		// 	constexpr resscr_t STD_DENOMINATOR = static_cast<resscr_t>( 50.000 );
		// 	constexpr resscr_t CUTOFF          = static_cast<resscr_t>(  0.001 );

		// 	const bool     is_suspicious  = ( ( arg_cond_evalue <= CUTOFF ) && ( arg_indp_evalue > CUTOFF ) && ( arg_cath_category != cath_id_score_category::DC_TYPE ) );
		// 	const bool     is_later_round =   arg_cath_category == cath_id_score_category::LATER_ROUND;
		// 	const resscr_t denom_mult     = is_suspicious  ? 4.0 :
		// 	                                is_later_round ? 2.0 :
		// 	                                                 1.0;
		// 	const resscr_t ratio          = ( arg_bitscore / ( STD_DENOMINATOR * denom_mult ) );
		// 	return ratio * ratio * ratio;
		// }

		void parse_domain_hits_table_file(read_and_process_mgr &,
		                                  const boost::filesystem::path &,
		                                  const bool &);

		void parse_domain_hits_table(read_and_process_mgr &,
		                             std::istream &,
		                             const bool &);

	} // namespace rslv
} // namespace cath

#endif
