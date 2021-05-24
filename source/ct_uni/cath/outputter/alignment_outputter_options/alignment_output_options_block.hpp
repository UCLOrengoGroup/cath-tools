/// \file
/// \brief The alignment_output_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS_ALIGNMENT_OUTPUT_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS_ALIGNMENT_OUTPUT_OPTIONS_BLOCK_HPP

#include <filesystem>
#include <iosfwd>

#include <boost/ptr_container/ptr_vector.hpp>

#include "cath/options/options_block/options_block.hpp"

// clang-format off
namespace cath { class display_spec; }
namespace cath::opts { class alignment_outputter; }
namespace cath::opts { class alignment_outputter_list; }
// clang-format on

namespace cath::opts {

	/// \brief TODOCUMENT
	class alignment_output_options_block final : public options_block {
	private:
		using alner = options_block;

		::std::filesystem::path aln_to_cath_aln_file;
		bool aln_to_cath_aln_stdout;
		::std::filesystem::path aln_to_fasta_file;
		bool aln_to_fasta_stdout;
		::std::filesystem::path aln_to_ssap_file;
		bool aln_to_ssap_stdout;
		::std::filesystem::path aln_to_html_file;
		bool aln_to_html_stdout;

		[[nodiscard]] std::unique_ptr<options_block> do_clone() const final;
		[[nodiscard]] std::string                    do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

		[[nodiscard]] ::std::filesystem::path get_aln_to_cath_aln_file() const;
		[[nodiscard]] bool                    get_aln_to_cath_aln_stdout() const;
		[[nodiscard]] ::std::filesystem::path get_aln_to_fasta_file() const;
		[[nodiscard]] bool                    get_aln_to_fasta_stdout() const;
		[[nodiscard]] ::std::filesystem::path get_aln_to_ssap_file() const;
		[[nodiscard]] bool                    get_aln_to_ssap_stdout() const;
		[[nodiscard]] ::std::filesystem::path get_aln_to_html_file() const;
		[[nodiscard]] bool                    get_aln_to_html_stdout() const;

	  public:
		[[nodiscard]] alignment_outputter_list get_alignment_outputters( const display_spec & ) const;

		[[nodiscard]] bool outputs_to_stdout() const;
	};

} // namespace cath::opts

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS_ALIGNMENT_OUTPUT_OPTIONS_BLOCK_HPP
