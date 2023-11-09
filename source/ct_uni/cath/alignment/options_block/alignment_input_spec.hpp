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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_OPTIONS_BLOCK_ALIGNMENT_INPUT_SPEC_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_OPTIONS_BLOCK_ALIGNMENT_INPUT_SPEC_HPP

#include <filesystem>

#include "cath/acquirer/alignment_acquirer/align_refining.hpp"
#include "cath/common/path_type_aliases.hpp"

namespace cath::opts {

	/// \brief Represent a specification for how alignments should be read in
	class alignment_input_spec final {
	private:
		/// \brief Whether to align based on matching residue names
		///
		/// This can be useful when aligning (consistently numbered) models of the same protein
		bool residue_name_align = DEFAULT_RESIDUE_NAME_ALIGN;
		
		/// \brief A file from which to read a FASTA alignment
		::std::filesystem::path fasta_alignment_file;

		/// \brief A file from which to read a legacy-SSAP-format alignment
		::std::filesystem::path ssap_alignment_file;

		/// \brief A file from which to read a CORA alignment
		::std::filesystem::path cora_alignment_file;

		/// \brief A file from which to read SSAP-scores format data to use to attempt to glue pairwise alignments together
		::std::filesystem::path ssap_scores_file;

		/// \brief A directory in which SSAPs should be performed and then their alignments glued together
		///        or (inner) nullopt for cath-tools to choose a directory to use
		///
		/// The outer optional<> is used to determine if the option has been specified
		/// The inner optional<> is used to determine if a specific directory should be used
		/// (rather than cath-tools choosing)
		path_opt_opt do_the_ssaps_dir;

		/// \brief How much refining should be done to the alignment
		align::align_refining refining = DEFAULT_REFINING;

	public:
		/// \brief The default value for whether to align based on matching residue names
		static constexpr bool DEFAULT_RESIDUE_NAME_ALIGN = false;

		/// \brief The default value for how much refining should be done to the alignment
		static constexpr align::align_refining DEFAULT_REFINING = align::align_refining::NO;

		alignment_input_spec() = default;
		explicit alignment_input_spec(const align::align_refining &);

		[[nodiscard]] const bool &                   get_residue_name_align() const;
		[[nodiscard]] const ::std::filesystem::path &get_fasta_alignment_file() const;
		[[nodiscard]] const ::std::filesystem::path &get_ssap_alignment_file() const;
		[[nodiscard]] const ::std::filesystem::path &get_cora_alignment_file() const;
		[[nodiscard]] const ::std::filesystem::path &get_ssap_scores_file() const;
		[[nodiscard]] const path_opt_opt &           get_do_the_ssaps_dir() const;
		[[nodiscard]] const align::align_refining &  get_refining() const;

		alignment_input_spec & set_residue_name_align(const bool &);
		alignment_input_spec & set_fasta_alignment_file(const ::std::filesystem::path &);
		alignment_input_spec & set_ssap_alignment_file(const ::std::filesystem::path &);
		alignment_input_spec & set_cora_alignment_file(const ::std::filesystem::path &);
		alignment_input_spec & set_ssap_scores_file(const ::std::filesystem::path &);
		alignment_input_spec & set_do_the_ssaps_dir(const path_opt &);
		alignment_input_spec & set_refining(const align::align_refining &);
	};

	size_t get_num_acquirers(const alignment_input_spec &);

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_OPTIONS_BLOCK_ALIGNMENT_INPUT_SPEC_HPP
