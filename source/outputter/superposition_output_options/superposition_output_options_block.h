/// \file
/// \brief The superposition_output_options_block class header

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

#ifndef SUPERPOSITION_OUTPUT_OPTIONS_BLOCK_H_INCLUDED
#define SUPERPOSITION_OUTPUT_OPTIONS_BLOCK_H_INCLUDED

#include <boost/ptr_container/ptr_vector.hpp>

#include "options/options_block/options_block.h"

#include <iosfwd>

namespace cath { namespace opts { class superposition_outputter; } }
namespace cath { namespace opts { class superposition_outputter_list; } }
namespace cath { class display_spec; }

namespace superposition_output_options_block_test_suite { struct parses_option_for_to_json_file; }
namespace superposition_output_options_block_test_suite { struct unparsed_has_no_json_file; }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class superposition_output_options_block final : public options_block {
		private:
			friend struct superposition_output_options_block_test_suite::parses_option_for_to_json_file;
			friend struct superposition_output_options_block_test_suite::unparsed_has_no_json_file;

			using super = options_block;

			static const std::string PO_SUP_FILE;
			static const std::string PO_SUP_FILES_DIR;
			static const std::string PO_SUP_TO_STDOUT;
			static const std::string PO_SUP_TO_PYMOL;
			static const std::string PO_PYMOL_PROGRAM;
			static const std::string PO_SUP_TO_PYMOL_FILE;
			static const std::string PO_SUP_TO_JSON_FILE;

			static const std::string DEFAULT_PYMOL_PROGRAM;

			boost::filesystem::path sup_to_pdb_file;
			boost::filesystem::path sup_to_pdb_files_dir;
			bool sup_to_stdout;
			bool sup_to_pymol;
			boost::filesystem::path sup_to_pymol_file;
			boost::filesystem::path pymol_program;
			boost::filesystem::path json_file;

			virtual std::unique_ptr<options_block> do_clone() const override final;
			virtual std::string do_get_block_name() const override final;
			virtual void do_add_visible_options_to_description(boost::program_options::options_description &) override final;
			virtual str_opt do_invalid_string(const boost::program_options::variables_map &) const override final;

			boost::filesystem::path get_sup_to_pdb_file() const;
			boost::filesystem::path get_sup_to_pdb_files_dir() const;
			bool get_sup_to_stdout() const;
			bool get_sup_to_pymol() const;
			boost::filesystem::path get_pymol_program() const;
			boost::filesystem::path get_sup_to_pymol_file() const;
			boost::filesystem::path get_json_file() const;

		public:
			virtual ~superposition_output_options_block() noexcept = default;

			/// \todo Consider adding a sister get_superposition_outputters() for getting non-display outputters
			///       without having to specify a display_spec
			superposition_outputter_list get_superposition_outputters(const display_spec &) const;
			bool outputs_to_stdout() const;
		};
	}
}

#endif
