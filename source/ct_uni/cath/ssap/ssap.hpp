/// \file
/// \brief Declarations of core SSAP functions

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 1989, Orengo Group, University College London
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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_SSAP_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_SSAP_HPP

#include <filesystem>
#include <iostream>
#include <optional>
#include <string>

#include "cath/alignment/align_type_aliases.hpp"
#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/common/path_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/ssap/compare_upper_cell_result.hpp"

// clang-format off
namespace cath { class entry_querier; }
namespace cath { class protein; }
namespace cath { class protein_source_file_set; }
namespace cath { class residue; }
namespace cath { class sec_struc; }
namespace cath { class selected_pair; }
namespace cath { class ssap_scores; }
namespace cath { struct clique; }
namespace cath::geom { class coord; }
namespace cath::opts { class cath_ssap_options; }
namespace cath::opts { class data_dirs_spec; }
namespace cath::opts { class old_ssap_options_block; }
// clang-format on

namespace cath {
	void reset_ssap_global_variables();

	void temp_set_global_run_counter(const ptrdiff_t &);

	ptrdiff_t temp_get_global_run_counter();

	prot_prot_pair read_protein_pair(const opts::cath_ssap_options &,
	                                 std::ostream & = std::cerr);

	prot_prot_pair read_protein_pair(const std::string &,
	                                 const chop::domain_opt &,
	                                 const std::string &,
	                                 const chop::domain_opt &,
	                                 const opts::data_dirs_spec &,
	                                 const protein_source_file_set &,
	                                 const path_opt &,
	                                 std::ostream & = std::cerr);

	void run_ssap(const opts::cath_ssap_options &,
	              std::ostream & = std::cout,
	              std::ostream & = std::cerr,
	              const ostream_ref_opt & = ::std::nullopt);

	void align_proteins(const protein &,
	                    const protein &,
	                    const opts::old_ssap_options_block &,
	                    const opts::data_dirs_spec &);

	ssap_scores fast_ssap(const protein &,
	                      const protein &,
	                      const opts::old_ssap_options_block &,
	                      const opts::data_dirs_spec &);

	std::pair<ssap_scores, align::alignment> compare(const protein &,
	                                                 const protein &,
	                                                 const size_t &,
	                                                 const entry_querier &,
	                                                 const opts::old_ssap_options_block &,
	                                                 const opts::data_dirs_spec &,
	                                                 const align::alignment_opt &);

	protein read_protein_data_from_ssap_options_files(const opts::data_dirs_spec &,
	                                                  const std::string &,
	                                                  const protein_source_file_set &,
	                                                  const path_opt &,
	                                                  const chop::domain_opt &,
	                                                  std::ostream & = std::cerr);

	clique read_clique_file(const ::std::filesystem::path &);

	void set_mask_matrix(const protein &,
	                     const protein &,
	                     const align::alignment_opt &,
	                     const path_opt &);

	void select_pairs(const protein &,
	                  const protein &,
	                  const size_t &,
	                  const entry_querier &);

	void update_best_pair_selections(std::deque<selected_pair> &,
	                                 const selected_pair &,
	                                 const size_t &);

	bool residues_have_similar_area_angle_props(const residue &,
	                                            const residue &);

	void populate_upper_score_matrix(const protein &,
	                                 const protein &,
	                                 const entry_querier &,
	                                 const bool &);

	compare_upper_cell_result compare_upper_cell(const protein &,
	                                             const protein &,
	                                             const size_t &,
	                                             const size_t &,
	                                             const entry_querier &,
	                                             const double &);

	score_type context_sec(const protein &,
	                       const protein &,
	                       const size_t &,
	                       const size_t &,
	                       const size_t &,
	                       const size_t &);


	ssap_scores calculate_log_score(const align::alignment &,
	                                const protein &,
	                                const protein &,
	                                const entry_querier &);

	double calculate_sequence_identity(const align::alignment &,
	                                   const protein &,
	                                   const protein &);

	bool save_ssap_scores(const align::alignment &,
	                      const protein &,
	                      const protein &,
	                      const ssap_scores &,
	                      const opts::old_ssap_options_block &,
	                      const opts::data_dirs_spec &);

	void save_zero_scores(const protein &,
	                      const protein &,
	                      const ptrdiff_t &);

	void print_ssap_scores(std::ostream &,
	                       const double &,
	                       const double &,
	                       const std::string &,
	                       const std::string &,
	                       const ptrdiff_t &,
	                       const bool &);

	size_doub_pair superpose(const protein &,
	                         const protein &,
	                         const align::alignment &,
	                         const opts::old_ssap_options_block &,
	                         const opts::data_dirs_spec &,
	                         const bool &);

	ssap_scores plot_aln(const protein &,
	                     const protein &,
	                     const size_t &,
	                     const entry_querier &,
	                     const align::alignment &,
	                     const opts::old_ssap_options_block &,
	                     const opts::data_dirs_spec &);

	::std::filesystem::path look_for_file_globally_then_locally(const std::string &,
	                                                            const std::string &,
	                                                            const ::std::filesystem::path &,
	                                                            const ::std::filesystem::path & = ::std::filesystem::current_path(),
	                                                            const bool & = false);

} // namespace cath

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_SSAP_HPP
