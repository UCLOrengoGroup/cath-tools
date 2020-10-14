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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_OPTIONS_BLOCK_ALIGNMENT_INPUT_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_OPTIONS_BLOCK_ALIGNMENT_INPUT_OPTIONS_BLOCK_HPP

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

			std::unique_ptr<options_block> do_clone() const final;
			std::string do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &,
			                                           const size_t &) final;
			str_opt do_invalid_string(const boost::program_options::variables_map &) const final;
			str_vec do_get_all_options_names() const final;

		public:
			static const std::string PO_RES_NAME_ALIGN;
			static const std::string PO_FASTA_ALIGN_INFILE;
			static const std::string PO_SSAP_ALIGN_INFILE;
			static const std::string PO_CORA_ALIGN_INFILE;
			static const std::string PO_SSAP_SCORE_INFILE;
			static const std::string PO_DO_THE_SSAPS;
			static const std::string PO_REFINING;

			alignment_input_options_block() = default;
			explicit alignment_input_options_block(const align::align_refining &);

			const alignment_input_spec & get_alignment_input_spec() const;
		};

		size_t get_num_acquirers(const alignment_input_options_block &);

	} // namespace opts
} // namespace cath

#endif
