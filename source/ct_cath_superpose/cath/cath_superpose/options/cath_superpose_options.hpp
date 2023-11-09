/// \file
/// \brief The cath_superpose_options class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CATH_SUPERPOSE_CATH_CATH_SUPERPOSE_OPTIONS_CATH_SUPERPOSE_OPTIONS_HPP
#define CATH_TOOLS_SOURCE_CT_CATH_SUPERPOSE_CATH_CATH_SUPERPOSE_OPTIONS_CATH_SUPERPOSE_OPTIONS_HPP

#include <iosfwd>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "cath/acquirer/alignment_acquirer/align_refining.hpp"
#include "cath/alignment/options_block/alignment_input_options_block.hpp"
#include "cath/alignment/options_block/alignment_input_spec.hpp"
#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/common/path_type_aliases.hpp"
#include "cath/display/options/display_options_block.hpp"
#include "cath/options/executable/executable_options.hpp"
#include "cath/options/options_block/ids_options_block.hpp"
#include "cath/options/options_block/pdb_input_options_block.hpp"
#include "cath/options/options_block/pdb_input_spec.hpp"
#include "cath/options/options_block/superposition_input_options_block.hpp"
#include "cath/outputter/alignment_outputter/alignment_outputter_list.hpp"
#include "cath/outputter/alignment_outputter_options/alignment_output_options_block.hpp"
#include "cath/outputter/superposition_output_options/superposition_output_options_block.hpp"
#include "cath/outputter/superposition_outputter/superposition_outputter_list.hpp"
#include "cath/superposition/options/align_regions_options_block.hpp"
#include "cath/superposition/options/superposition_content_options_block.hpp"

// clang-format off
namespace cath::align { class alignment; }
namespace cath::align { class alignment_acquirer; }
namespace cath::file { class strucs_context; }
namespace cath::opts { class pdbs_acquirer; }
namespace cath::opts { class selection_policy_acquirer; }
namespace cath::sup { class superposition_context; }
// clang-format on

namespace cath::opts {

	/// \brief TODOCUMENT
	///
	/// \todo Add options to allow identifiable PDBs to be loaded from a directory
	/// \todo Sort out the interaction between loading PDBs and loading an alignment
	/// \todo Add options to write/read a superposition
	///
	/// \todo Abstract out duplication between cath_refine_align_options / cath_superpose_options
	class cath_superpose_options final : public executable_options {
	private:
		using super = executable_options;

		/// \brief The options_block for the alignment input options
		alignment_input_options_block            the_alignment_input_ob{ align::align_refining::NO };

		/// \brief The options_block for the superposition input options
		superposition_input_options_block        the_superposition_input_ob;

		/// \brief The options_block for the IDs options
		ids_options_block                        the_ids_ob;

		/// \brief The options_block for the pdb input options
		pdb_input_options_block                  the_pdb_input_ob;

		/// \brief The align_regions_options_block for align regions options
		align_regions_options_block              the_align_regions_ob;

		/// \brief The options_block for the alignment output options
		alignment_output_options_block           the_alignment_output_ob;

		/// \brief The options_block for the superposition output options
		superposition_output_options_block       the_superposition_output_ob;

		/// \brief The options_block for the display options
		display_options_block                    the_display_ob;

		/// \brief The specification of what should be included in the superposition
		sup::superposition_content_options_block the_content_ob;

		[[nodiscard]] std::string_view do_get_program_name() const final;
		[[nodiscard]] str_opt          do_get_error_or_help_string() const final;

		[[nodiscard]] std::string do_get_help_prefix_string() const final;
		[[nodiscard]] std::string do_get_help_suffix_string() const final;
		[[nodiscard]] std::string do_get_overview_string() const final;

	  public:
		cath_superpose_options();

		void check_ok_to_use() const;

		[[nodiscard]] const path_opt &get_json_sup_infile() const;

		[[nodiscard]] selection_policy_acquirer get_selection_policy_acquirer() const;

		[[nodiscard]] const alignment_input_spec &get_alignment_input_spec() const;
		[[nodiscard]] const str_vec &             get_ids() const;
		[[nodiscard]] const pdb_input_spec &      get_pdb_input_spec() const;

		[[nodiscard]] alignment_outputter_list     get_alignment_outputters() const;
		[[nodiscard]] superposition_outputter_list get_superposition_outputters( const default_supn_outputter & ) const;

		[[nodiscard]] const chop::domain_vec &get_domains() const;

		/// \brief The name of this program
		static constexpr ::std::string_view PROGRAM_NAME{ "cath-superpose" };
	};

	std::unique_ptr<const align::alignment_acquirer> get_alignment_acquirer(const cath_superpose_options &);
	align::align_refining get_align_refining(const cath_superpose_options &);
	std::pair<align::alignment, size_size_pair_vec> get_alignment_and_spanning_tree(const cath_superpose_options &,
	                                                                                const file::strucs_context &,
	                                                                                const align::align_refining &);

	std::unique_ptr<const pdbs_acquirer> get_pdbs_acquirer(const cath_superpose_options &);
	file::strucs_context get_pdbs_and_names(const cath_superpose_options &,
	                                        std::istream &,
	                                        const bool &);

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_CATH_SUPERPOSE_CATH_CATH_SUPERPOSE_OPTIONS_CATH_SUPERPOSE_OPTIONS_HPP
