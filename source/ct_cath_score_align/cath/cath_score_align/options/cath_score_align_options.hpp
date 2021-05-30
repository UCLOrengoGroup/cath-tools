/// \file
/// \brief The cath_score_align_options class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CATH_SCORE_ALIGN_CATH_CATH_SCORE_ALIGN_OPTIONS_CATH_SCORE_ALIGN_OPTIONS_HPP
#define _CATH_TOOLS_SOURCE_CT_CATH_SCORE_ALIGN_CATH_CATH_SCORE_ALIGN_OPTIONS_CATH_SCORE_ALIGN_OPTIONS_HPP

#include <iosfwd>
#include <string_view>
#include <vector>

#include "cath/alignment/options_block/alignment_input_options_block.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/options/executable/executable_options.hpp"
#include "cath/options/options_block/pdb_input_options_block.hpp"

// clang-format off
namespace cath::align { class alignment_acquirer; }
namespace cath::file { class strucs_context; }
namespace cath::opts { class pdbs_acquirer; }
// clang-format on

namespace cath::opts {

	/// \brief TODOCUMENT
	class cath_score_align_options final : public executable_options {
	private:
		using super = executable_options;

		/// \brief TODOCUMENT
		alignment_input_options_block the_alignment_input_options_block{ align::align_refining::NO };

		/// \brief TODOCUMENT
		pdb_input_options_block       the_pdb_input_options_block;

		[[nodiscard]] std::string_view do_get_program_name() const final;
		[[nodiscard]] str_opt          do_get_error_or_help_string() const final;

		[[nodiscard]] std::string do_get_help_prefix_string() const final;
		[[nodiscard]] std::string do_get_help_suffix_string() const final;
		[[nodiscard]] std::string do_get_overview_string() const final;

		void check_ok_to_use() const;

	public:
		cath_score_align_options();

		[[nodiscard]] const pdb_input_spec &      get_pdb_input_spec() const;
		[[nodiscard]] const alignment_input_spec &get_alignment_input_spec() const;

		/// \brief The name of the program that uses this executable_options
		static constexpr ::std::string_view PROGRAM_NAME{ "cath-score-align" };
	};

	std::unique_ptr<const align::alignment_acquirer> get_alignment_acquirer(const cath_score_align_options &);
	std::unique_ptr<const pdbs_acquirer> get_pdbs_acquirer(const cath_score_align_options &);
	file::strucs_context get_pdbs_and_names(const cath_score_align_options &,
	                                        std::istream &,
	                                        const bool &);

} // namespace cath::opts

#endif // _CATH_TOOLS_SOURCE_CT_CATH_SCORE_ALIGN_CATH_CATH_SCORE_ALIGN_OPTIONS_CATH_SCORE_ALIGN_OPTIONS_HPP
