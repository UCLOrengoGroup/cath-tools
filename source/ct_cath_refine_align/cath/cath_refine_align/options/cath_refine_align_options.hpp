/// \file
/// \brief The cath_refine_align_options class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CATH_REFINE_ALIGN_CATH_CATH_REFINE_ALIGN_OPTIONS_CATH_REFINE_ALIGN_OPTIONS_HPP
#define _CATH_TOOLS_SOURCE_CT_CATH_REFINE_ALIGN_CATH_CATH_REFINE_ALIGN_OPTIONS_CATH_REFINE_ALIGN_OPTIONS_HPP

#include "cath/alignment/options_block/alignment_input_options_block.hpp"
#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/display/options/display_options_block.hpp"
#include "cath/options/executable/executable_options.hpp"
#include "cath/options/options_block/ids_options_block.hpp"
#include "cath/options/options_block/pdb_input_options_block.hpp"
#include "cath/outputter/alignment_outputter_options/alignment_output_options_block.hpp"
#include "cath/outputter/superposition_output_options/superposition_output_options_block.hpp"
#include "cath/superposition/options/align_regions_options_block.hpp"

#include <iosfwd>
#include <vector>

namespace cath { namespace align { class alignment_acquirer; } }
namespace cath { namespace file { class strucs_context; } }
namespace cath { namespace opts { class pdbs_acquirer; } }
namespace cath { namespace opts { class selection_policy_acquirer; } }
namespace cath { namespace opts { class superposition_outputter; } }
namespace cath { namespace sup { class superposition_context; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		///
		/// \todo Abstract out duplication between cath_refine_align_options / cath_superpose_options
		class cath_refine_align_options final : public executable_options {
		private:
			using super = executable_options;

			static const std::string STANDARD_USAGE_ERROR_STRING;

			/// \brief TODOCUMENT
			alignment_input_options_block      the_alignment_input_options_block{ align::align_refining::HEAVY };

			/// \brief The options_block for the IDs options
			ids_options_block                  the_ids_ob;

			/// \brief TODOCUMENT
			pdb_input_options_block            the_pdb_input_options_block;

			/// \brief The align_regions_options_block for align regions options
			align_regions_options_block        the_align_regions_ob;

			/// \brief TODOCUMENT
			alignment_output_options_block     the_alignment_output_options_block;

			/// \brief TODOCUMENT
			superposition_output_options_block the_superposition_output_options_block;

			/// \brief TODOCUMENT
			display_options_block              the_display_options_block;

			/// \brief The specification of what should be included in the superposition
			sup::superposition_content_spec    the_content_spec;

			[[nodiscard]] std::string do_get_program_name() const final;
			[[nodiscard]] str_opt     do_get_error_or_help_string() const final;

			[[nodiscard]] std::string do_get_help_prefix_string() const final;
			[[nodiscard]] std::string do_get_help_suffix_string() const final;
			[[nodiscard]] std::string do_get_overview_string() const final;

			void check_ok_to_use() const;

		public:
			cath_refine_align_options();

			[[nodiscard]] const str_vec &             get_ids() const;
			[[nodiscard]] const pdb_input_spec &      get_pdb_input_spec() const;
			[[nodiscard]] const alignment_input_spec &get_alignment_input_spec() const;
			[[nodiscard]] alignment_outputter_list    get_alignment_outputters() const;
			[[nodiscard]] superposition_outputter_list get_superposition_outputters( const default_supn_outputter & ) const;

			[[nodiscard]] const chop::domain_vec &get_domains() const;

			static const std::string PROGRAM_NAME;
		};

		std::unique_ptr<const align::alignment_acquirer> get_alignment_acquirer(const cath_refine_align_options &);
		align::align_refining get_align_refining(const cath_refine_align_options &);
		std::unique_ptr<const pdbs_acquirer> get_pdbs_acquirer(const cath_refine_align_options &);
		selection_policy_acquirer get_selection_policy_acquirer(const cath_refine_align_options &);

		file::strucs_context get_pdbs_and_names(const cath_refine_align_options &,
		                                        std::istream &,
		                                        const bool &);
	} // namespace opts
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CATH_REFINE_ALIGN_CATH_CATH_REFINE_ALIGN_OPTIONS_CATH_REFINE_ALIGN_OPTIONS_HPP
