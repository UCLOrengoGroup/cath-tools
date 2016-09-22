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

#ifndef CATH_SUPERPOSE_OPTIONS_H_INCLUDED
#define CATH_SUPERPOSE_OPTIONS_H_INCLUDED

#include "options/executable/executable_options.h"
#include "options/options_block/alignment_input_options_block.h"
#include "options/options_block/alignment_output_options_block.h"
#include "options/options_block/display_options_block.h"
#include "options/options_block/ids_options_block.h"
#include "options/options_block/misc_help_version_options_block.h"
#include "options/options_block/pdb_input_options_block.h"
#include "options/options_block/superposition_output_options_block.h"
#include "options/outputter/alignment_outputter/alignment_outputter_list.h"
#include "options/outputter/superposition_outputter/superposition_outputter_list.h"

#include <iosfwd>
#include <memory>
#include <string>

namespace boost { namespace program_options { class options_description; } }

namespace cath { namespace opts { class pdbs_acquirer; } }
namespace cath { namespace opts { class selection_policy_acquirer; } }
namespace cath { namespace sup { class superposition_context; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		///
		/// \todo Add options to allow identifiable PDBs to be loaded from a directory
		/// \todo Sort out the interaction between loading PDBs and loading an alignment
		/// \todo Add options to write/read a superposition
		class cath_superpose_options final : public executable_options {
		private:
			using super = executable_options;

			static const std::string STANDARD_USAGE_ERROR_STRING;

			misc_help_version_options_block    the_misc_options_block;
			alignment_input_options_block      the_alignment_input_options_block;
			ids_options_block                  the_ids_options_block;
			pdb_input_options_block            the_pdb_input_options_block;
			alignment_output_options_block     the_alignment_output_options_block;
			superposition_output_options_block the_superposition_output_options_block;
			display_options_block              the_display_options_block;

			virtual std::string do_get_program_name() const override final;
			virtual std::string do_update_error_or_help_string(const boost::program_options::options_description &) const override final;

			static std::string get_help_prefix_string();
			static std::string get_help_suffix_string();

			void check_ok_to_use() const;

			std::unique_ptr<const pdbs_acquirer> get_pdbs_acquirer() const;
			selection_policy_acquirer get_selection_policy_acquirer() const;

		public:
			cath_superpose_options();
			virtual ~cath_superpose_options() noexcept = default;

			sup::superposition_context get_superposition_context(std::istream &,
			                                                     std::ostream &) const;
			alignment_outputter_list get_alignment_outputters() const;
			superposition_outputter_list get_superposition_outputters() const;

			static std::string get_overview_string();

			static const std::string PROGRAM_NAME;
		};

	}
}

#endif
