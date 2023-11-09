/// \file
/// \brief The ssap_and_prc class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_SSAP_AND_PRC_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_SSAP_AND_PRC_HPP

#include <optional>

#include "cath/file/prc_scores_file/prc_scores_entry.hpp"
#include "cath/file/ssap_scores_file/ssap_scores_entry.hpp"

// clang-format off
namespace cath::file { class prc_scores_entry; }
namespace cath::file { class ssap_scores_entry; }
namespace cath::score { class rbf_model; }
// clang-format on

namespace cath::homcheck {

	/// \brief Represent a SSAP result and a PRC result for the same query/match pair
	class ssap_and_prc final {
	private:
		/// \brief The SSAP result
		file::ssap_scores_entry the_ssap;

		/// \brief The PRC result
		file::prc_scores_entry the_prc;

		/// \brief The magic function score calculated from the SSAP and PRC results
		///
		/// This gets populated by the ctor
		double magic_function_score;

		/// \brief The SVM score for this SSAP and PRC pair which may or may not get populated
		::std::optional<double> svm_score;

	public:
		ssap_and_prc(file::ssap_scores_entry,
		             file::prc_scores_entry);

		void calculate_svm_score(const score::rbf_model &);

		[[nodiscard]] const std::string &get_query_id() const;
		[[nodiscard]] const std::string &get_match_id() const;

		[[nodiscard]] const file::ssap_scores_entry &get_ssap() const;
		[[nodiscard]] const file::prc_scores_entry & get_prc() const;
		[[nodiscard]] const double &                 get_magic_function_score() const;

		[[nodiscard]] const ::std::optional<double> &get_svm_score_opt() const;
	};

	double magic_function(const file::ssap_scores_entry &,
	                      const file::prc_scores_entry &);

	bool has_svm_score(const ssap_and_prc &);
	const double & get_svm_score(const ssap_and_prc &);

	const size_t & get_ssap_length_1(const ssap_and_prc &);
	const size_t & get_ssap_length_2(const ssap_and_prc &);
	const double & get_ssap_score(const ssap_and_prc &);
	const size_t & get_ssap_num_equivs(const ssap_and_prc &);
	const double & get_ssap_overlap_pc(const ssap_and_prc &);
	const double & get_ssap_seq_id_pc(const ssap_and_prc &);
	const double & get_ssap_rmsd(const ssap_and_prc &);
	const size_t & get_prc_start_1(const ssap_and_prc &);
	const size_t & get_prc_end_1(const ssap_and_prc &);
	const size_t & get_prc_length_1(const ssap_and_prc &);
	const size_t & get_prc_hit_num(const ssap_and_prc &);
	const size_t & get_prc_start_2(const ssap_and_prc &);
	const size_t & get_prc_end_2(const ssap_and_prc &);
	const size_t & get_prc_length_2(const ssap_and_prc &);
	const double & get_prc_simple(const ssap_and_prc &);
	const double & get_prc_reverse(const ssap_and_prc &);
	const double & get_prc_evalue(const ssap_and_prc &);

	std::string to_string(const ssap_and_prc &);

	std::ostream & operator<<(std::ostream &,
	                          const ssap_and_prc &);

} // namespace cath::homcheck

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_SSAP_AND_PRC_HPP
