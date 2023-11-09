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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS_SUPERPOSITION_OUTPUT_OPTIONS_BLOCK_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS_SUPERPOSITION_OUTPUT_OPTIONS_BLOCK_HPP

#include <filesystem>
#include <iosfwd>
#include <string_view>

#include <boost/ptr_container/ptr_vector.hpp>

#include "cath/options/options_block/options_block.hpp"
#include "cath/superposition/superposition_content_spec.hpp"

// clang-format off
namespace cath { class display_spec; }
namespace cath::opts { class superposition_outputter; }
namespace cath::opts { class superposition_outputter_list; }
namespace cath::sup { class superposition_content_spec; }
// clang-format on

// clang-format off
namespace superposition_output_options_block_test_suite { struct parses_option_for_to_json_file; }
namespace superposition_output_options_block_test_suite { struct unparsed_has_no_json_file; }
// clang-format on

namespace cath::opts {

	/// \brief What superposition outputter (if any) should be provided if none is explicitly specified
	enum class default_supn_outputter : bool {
		PYMOL, ///< The pymol_view_superposition_outputter
		NONE   ///< No outputter
	};

	/// \brief TODOCUMENT
	class superposition_output_options_block final : public options_block {
	private:
		friend struct ::superposition_output_options_block_test_suite::parses_option_for_to_json_file;
		friend struct ::superposition_output_options_block_test_suite::unparsed_has_no_json_file;

		using super = options_block;

		::std::filesystem::path sup_to_pdb_file;
		::std::filesystem::path sup_to_pdb_files_dir;
		bool sup_to_stdout;
		bool sup_to_pymol;
		::std::filesystem::path sup_to_pymol_file;
		::std::filesystem::path pymol_program;
		::std::filesystem::path json_file;

		[[nodiscard]] std::unique_ptr<options_block> do_clone() const final;
		[[nodiscard]] std::string                    do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

		[[nodiscard]] ::std::filesystem::path get_sup_to_pdb_file() const;
		[[nodiscard]] ::std::filesystem::path get_sup_to_pdb_files_dir() const;
		[[nodiscard]] bool                    get_sup_to_stdout() const;
		[[nodiscard]] bool                    get_sup_to_pymol() const;
		[[nodiscard]] ::std::filesystem::path get_pymol_program() const;
		[[nodiscard]] ::std::filesystem::path get_sup_to_pymol_file() const;
		[[nodiscard]] ::std::filesystem::path get_json_file() const;

	  public:
		/// \todo Consider adding a sister get_superposition_outputters() for getting outputters
		///       that don't require a display_spec / superposition_content_spec
		[[nodiscard]] superposition_outputter_list get_superposition_outputters(
		  const display_spec &,
		  const sup::superposition_content_spec & = sup::superposition_content_spec{},
		  const default_supn_outputter &          = default_supn_outputter::NONE ) const;
		[[nodiscard]] bool outputs_to_stdout() const;

		static constexpr ::std::string_view DEFAULT_PYMOL_PROGRAM{ "pymol" };
		static constexpr ::std::string_view PO_PYMOL_PROGRAM{ "pymol-program" };
		static constexpr ::std::string_view PO_SUP_FILE{ "sup-to-pdb-file" };
		static constexpr ::std::string_view PO_SUP_FILES_DIR{ "sup-to-pdb-files-dir" };
		static constexpr ::std::string_view PO_SUP_TO_JSON_FILE{ "sup-to-json-file" };
		static constexpr ::std::string_view PO_SUP_TO_PYMOL{ "sup-to-pymol" };
		static constexpr ::std::string_view PO_SUP_TO_PYMOL_FILE{ "sup-to-pymol-file" };
		static constexpr ::std::string_view PO_SUP_TO_STDOUT{ "sup-to-stdout" };
	};

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS_SUPERPOSITION_OUTPUT_OPTIONS_BLOCK_HPP
