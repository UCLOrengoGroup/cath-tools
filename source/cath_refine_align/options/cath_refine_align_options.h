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

#ifndef CATH_REFINE_ALIGN_OPTIONS_H_INCLUDED
#define CATH_REFINE_ALIGN_OPTIONS_H_INCLUDED

#include "common/type_aliases.h"
#include "display/options/display_options_block.h"
#include "options/executable/executable_options.h"
#include "options/options_block/alignment_input_options_block.h"
#include "options/options_block/pdb_input_options_block.h"
#include "outputter/alignment_outputter_options/alignment_output_options_block.h"
#include "outputter/superposition_output_options/superposition_output_options_block.h"

#include <iosfwd>
#include <vector>

namespace cath { namespace opts { class alignment_acquirer; } }
namespace cath { namespace opts { class pdbs_acquirer; } }
namespace cath { namespace opts { class selection_policy_acquirer; } }
namespace cath { namespace opts { class superposition_outputter; } }
namespace cath { namespace sup { class superposition_context; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class cath_refine_align_options final : public executable_options {
		private:
			using super = executable_options;

			static const std::string STANDARD_USAGE_ERROR_STRING;

			/// \brief TODOCUMENT
			alignment_input_options_block      the_alignment_input_options_block;

			/// \brief TODOCUMENT
			pdb_input_options_block            the_pdb_input_options_block;

			/// \brief TODOCUMENT
			alignment_output_options_block     the_alignment_output_options_block;

			/// \brief TODOCUMENT
			superposition_output_options_block the_superposition_output_options_block;

			/// \brief TODOCUMENT
			display_options_block              the_display_options_block;

			virtual std::string do_get_program_name() const override final;
			virtual str_opt do_get_error_or_help_string() const override final;

			virtual std::string do_get_help_prefix_string() const override final;
			virtual std::string do_get_help_suffix_string() const override final;
			virtual std::string do_get_overview_string() const override final;

			void check_ok_to_use() const;

		public:
			cath_refine_align_options();
			virtual ~cath_refine_align_options() noexcept = default;

			const pdb_input_spec & get_pdb_input_spec() const;
			const alignment_input_spec & get_alignment_input_spec() const;
			alignment_outputter_list get_alignment_outputters() const;
			superposition_outputter_list get_superposition_outputters() const;

			static const std::string PROGRAM_NAME;
		};

		std::unique_ptr<const alignment_acquirer> get_alignment_acquirer(const cath_refine_align_options &);
		std::unique_ptr<const pdbs_acquirer> get_pdbs_acquirer(const cath_refine_align_options &);
		selection_policy_acquirer get_selection_policy_acquirer(const cath_refine_align_options &);
	}
}

#endif
