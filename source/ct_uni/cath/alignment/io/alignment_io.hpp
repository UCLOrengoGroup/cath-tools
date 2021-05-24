/// \file
/// \brief The alignment io class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_IO_ALIGNMENT_IO_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_IO_ALIGNMENT_IO_HPP

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>

#include "cath/alignment/align_type_aliases.hpp"
#include "cath/alignment/alignment.hpp"
#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/structure/structure_type_aliases.hpp"

// clang-format off
namespace cath { class protein; }
namespace cath { class protein_list; }
namespace cath::file { class name_set_list; }
namespace cath::file { class pdb; }
namespace cath::file { class pdb_list; }
// clang-format on

namespace cath::align {

	/// \todo Organise these into three class hierarchies: alignment_format, alignment_reader and alignment_writer
	///       a la the chopping hierarchies
	alignment read_alignment_from_cath_ssap_legacy_format(const ::std::filesystem::path &,
	                                                      const protein &,
	                                                      const protein &,
	                                                      const ostream_ref_opt & = std::ref( std::cerr ) );

	alignment read_alignment_from_cath_ssap_legacy_format(std::istream &,
	                                                      const protein &,
	                                                      const protein &,
	                                                      const ostream_ref_opt & = std::ref( std::cerr ) );

	alignment read_alignment_from_cath_ssap_legacy_format(std::istream &,
	                                                      const file::pdb &,
	                                                      const file::pdb &,
	                                                      const ostream_ref_opt & = std::ref( std::cerr ) );

	alignment read_alignment_from_cath_ssap_legacy_format(std::istream &,
	                                                      const residue_id_vec &,
	                                                      const residue_id_vec &,
	                                                      const ostream_ref_opt & = std::ref( std::cerr ) );

	alignment read_alignment_from_cath_cora_legacy_format(std::istream &,
	                                                      const file::pdb_list &,
	                                                      const ostream_ref_opt & = std::ref( std::cerr ) );

	str_str_pair_vec read_ids_and_sequences_from_fasta(std::istream &);

	aln_posn_opt_vec align_sequence_to_amino_acids(const std::string &,
	                                               const amino_acid_vec &,
	                                               const std::string &,
	                                               std::ostream & = std::cerr);

	/// \todo Write code to make it easy to parse from proteins (and to optionally rescore whilst doing so)
	alignment read_alignment_from_fasta_file(const ::std::filesystem::path &,
	                                         const file::pdb_list &,
	                                         std::ostream & = std::cerr);

	alignment read_alignment_from_fasta_file(const ::std::filesystem::path &,
	                                         const file::pdb_list &,
	                                         const str_vec &,
	                                         std::ostream & = std::cerr);

	alignment read_alignment_from_fasta(std::istream &,
	                                    const file::pdb_list &,
	                                    std::ostream & = std::cerr);

	alignment read_alignment_from_fasta_file(const ::std::filesystem::path &,
	                                         const protein_list &,
	                                         std::ostream & = std::cerr);

	alignment read_alignment_from_fasta_file(const ::std::filesystem::path &,
	                                         const protein_list &,
	                                         const str_vec &,
	                                         std::ostream & = std::cerr);

	alignment read_alignment_from_fasta(std::istream &,
	                                    const protein_list &,
	                                    std::ostream & = std::cerr);

	alignment read_alignment_from_fasta(std::istream &,
	                                    const amino_acid_vec_vec &,
	                                    const str_vec &,
	                                    std::ostream & = std::cerr);

	aln_posn_opt search_for_residue_in_residue_ids(const size_t &,
	                                               const residue_id_vec &,
	                                               const char &,
	                                               const residue_name &,
	                                               const ostream_ref_opt &);

	void write_alignment_as_cath_ssap_legacy_format(const ::std::filesystem::path &,
	                                                const alignment &,
	                                                const protein &,
	                                                const protein &,
	                                                const chop::region_vec_opt & = ::std::nullopt,
	                                                const chop::region_vec_opt & = ::std::nullopt);

	std::ostream & output_alignment_to_cath_ssap_legacy_format(std::ostream &,
	                                                           const alignment &,
	                                                           const protein &,
	                                                           const protein &,
	                                                           const chop::region_vec_opt & = ::std::nullopt,
	                                                           const chop::region_vec_opt & = ::std::nullopt);

	std::string to_cath_ssap_legacy_format_alignment_string(const alignment &,
	                                                        const protein &,
	                                                        const protein &,
	                                                        const chop::region_vec_opt & = ::std::nullopt,
	                                                        const chop::region_vec_opt & = ::std::nullopt);

	void write_alignment_as_fasta_alignment(const ::std::filesystem::path &,
	                                        const alignment &,
	                                        const protein_list &);

	std::ostream & write_alignment_as_fasta_alignment(std::ostream &,
	                                                  const alignment &,
	                                                  const protein_list &);

	std::ostream & write_alignment_as_fasta_alignment(std::ostream &,
	                                                  const alignment &,
	                                                  const file::pdb_list &,
	                                                  const file::name_set_list &);

	std::string alignment_as_fasta_string(const alignment &,
	                                      const protein_list &);

	std::string alignment_as_fasta_string(const alignment &,
	                                      const file::pdb_list &,
	                                      const file::name_set_list &);

	std::string alignment_as_fasta_string(const alignment &,
	                                      const file::pdb_list &,
	                                      const str_vec &);

} // namespace cath::align

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_IO_ALIGNMENT_IO_HPP
