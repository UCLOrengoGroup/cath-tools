/// \file
/// \brief The alignment_input_spec class header

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

#ifndef ALIGNMENT_INPUT_SPEC_H_INCLUDED
#define ALIGNMENT_INPUT_SPEC_H_INCLUDED

#include <boost/filesystem/path.hpp>

namespace cath {
	namespace opts {

		/// \brief Represent a specification for how alignments should be read in
		class alignment_input_spec final {
		private:
			/// \brief Whether to align based on matching residue names
			///
			/// This can be useful when aligning (consistently numbered) models of the same protein
			bool residue_name_align = DEFAULT_RESIDUE_NAME_ALIGN;
			
			/// \brief A file from which to read a FASTA alignment
			boost::filesystem::path fasta_alignment_file;

			/// \brief A file from which to read a legacy-SSAP-format alignment
			boost::filesystem::path ssap_alignment_file;

			/// \brief A file from which to read a CORA alignment
			boost::filesystem::path cora_alignment_file;

			/// \brief A file from which to read SSAP-scores format data to use to attempt to glue pairwise alignments together
			boost::filesystem::path ssap_scores_file;

		public:
			/// \brief The default value for whether to align based on matching residue names
			static constexpr bool DEFAULT_RESIDUE_NAME_ALIGN = false;

			const bool & get_residue_name_align() const;
			const boost::filesystem::path & get_fasta_alignment_file() const;
			const boost::filesystem::path & get_ssap_alignment_file() const;
			const boost::filesystem::path & get_cora_alignment_file() const;
			const boost::filesystem::path & get_ssap_scores_file() const;

			alignment_input_spec & set_residue_name_align(const bool &);
			alignment_input_spec & set_fasta_alignment_file(const boost::filesystem::path &);
			alignment_input_spec & set_ssap_alignment_file(const boost::filesystem::path &);
			alignment_input_spec & set_cora_alignment_file(const boost::filesystem::path &);
			alignment_input_spec & set_ssap_scores_file(const boost::filesystem::path &);
		};

		size_t get_num_acquirers(const alignment_input_spec &);

	}
}

#endif
