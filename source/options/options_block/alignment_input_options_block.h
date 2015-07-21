/// \file
/// \brief The alignment_input_options_block class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef ALIGNMENT_INPUT_OPTIONS_BLOCK_H_INCLUDED
#define ALIGNMENT_INPUT_OPTIONS_BLOCK_H_INCLUDED

#include <boost/ptr_container/ptr_vector.hpp>

#include "options/options_block/options_block.h"

namespace cath { namespace opts { class alignment_acquirer; } }
namespace cath { namespace file { class pdb_list; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class alignment_input_options_block final : public options_block {
		private:
			using super = options_block;

			static const std::string PO_RES_NAME_ALIGN;
			static const std::string PO_FASTA_ALIGN_INFILE;
			static const std::string PO_SSAP_ALIGN_INFILE;
			static const std::string PO_CORA_ALIGN_INFILE;
			static const std::string PO_SSAP_SCORE_INFILE;

			bool residue_name_align;
			boost::filesystem::path fasta_alignment_file;
			boost::filesystem::path ssap_alignment_file;
			boost::filesystem::path cora_alignment_file;
			boost::filesystem::path ssap_scores_file;

			virtual std::unique_ptr<options_block> do_clone() const override final;
			virtual std::string do_get_block_name() const override final;
			virtual void do_add_visible_options_to_description(boost::program_options::options_description &) override final;
			virtual opt_str do_invalid_string() const override final;

			boost::filesystem::path get_fasta_alignment_file() const;
			boost::filesystem::path get_ssap_alignment_file() const;
			boost::filesystem::path get_cora_alignment_file() const;

		public:
			virtual ~alignment_input_options_block() noexcept = default;

			bool get_residue_name_align() const;
			boost::filesystem::path get_ssap_scores_file() const;

			boost::ptr_vector<alignment_acquirer> get_alignment_acquirers() const;
		};

	}
}

#endif
