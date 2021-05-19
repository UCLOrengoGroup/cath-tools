/// \file
/// \brief The alignment_input_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_OPTIONS_BLOCK_ALIGNMENT_INPUT_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_OPTIONS_BLOCK_ALIGNMENT_INPUT_OPTIONS_BLOCK_HPP

#include <string_view>

#include "cath/alignment/options_block/alignment_input_spec.hpp"
#include "cath/options/options_block/options_block.hpp"

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class alignment_input_options_block final : public options_block {
		private:
			using super = options_block;

			/// \brief The alignment_input_spec to be configured by this options block
			alignment_input_spec the_alignment_input_spec;

			[[nodiscard]] std::unique_ptr<options_block> do_clone() const final;
			[[nodiscard]] std::string                    do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &,
			                                           const size_t &) final;
			[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
			[[nodiscard]] str_view_vec do_get_all_options_names() const final;

		  public:
			alignment_input_options_block() = default;
			explicit alignment_input_options_block(const align::align_refining &);

			[[nodiscard]] const alignment_input_spec &get_alignment_input_spec() const;

			/// \brief The option name for whether to align based on matching residue names
			static constexpr ::std::string_view PO_RES_NAME_ALIGN{ "res-name-align" };

			/// \brief The option name for a file from which to read a FASTA alignment
			static constexpr ::std::string_view PO_FASTA_ALIGN_INFILE{ "fasta-aln-infile" };

			/// \brief The option name for a file from which to read a legacy-SSAP-format alignment
			static constexpr ::std::string_view PO_SSAP_ALIGN_INFILE{ "ssap-aln-infile" };

			/// \brief The option name for a file from which to read a CORA alignment
			static constexpr ::std::string_view PO_CORA_ALIGN_INFILE{ "cora-aln-infile" };

			/// \brief The option name for a file from which to read SSAP-scores format data to use to attempt to glue pairwise alignments together
			static constexpr ::std::string_view PO_SSAP_SCORE_INFILE{ "ssap-scores-infile" };

			/// \brief The option name for a directory in which to do the necessary SSAPs and then use the scores to glue the resulting alignments together
			static constexpr ::std::string_view PO_DO_THE_SSAPS{ "do-the-ssaps" };

			/// \brief The option name for how much refining should be done to the alignment
			static constexpr ::std::string_view PO_REFINING{ "align-refining" };
		};

		size_t get_num_acquirers(const alignment_input_options_block &);

	} // namespace opts
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_OPTIONS_BLOCK_ALIGNMENT_INPUT_OPTIONS_BLOCK_HPP
