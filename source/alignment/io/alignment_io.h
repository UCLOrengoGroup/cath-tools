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

#ifndef _CATH_TOOLS_SOURCE_ALIGNMENT_IO_ALIGNMENT_IO_H
#define _CATH_TOOLS_SOURCE_ALIGNMENT_IO_ALIGNMENT_IO_H

#include <boost/filesystem/path.hpp>

#include "alignment/align_type_aliases.h"
#include "alignment/alignment.h"
#include "common/type_aliases.h"
#include "structure/structure_type_aliases.h"

#include <cstddef>
#include <iostream>
#include <string>

namespace cath { namespace file { class pdb; } }
namespace cath { namespace file { class pdb_list; } }
namespace cath { class protein; }
namespace cath { class protein_list; }


namespace cath {

	namespace align {

		/// \todo Organise these into three class hierarchies: alignment_format, alignment_reader and alignment_writer
		///       a la the chopping hierarchies

		alignment read_alignment_from_cath_ssap_legacy_format(const boost::filesystem::path &,
		                                                      const protein &,
		                                                      const protein &,
		                                                      std::ostream &arg_stderr = std::cerr);

		alignment read_alignment_from_cath_ssap_legacy_format(std::istream &,
		                                                      const protein &,
		                                                      const protein &,
		                                                      std::ostream &arg_stderr = std::cerr);

		alignment read_alignment_from_cath_ssap_legacy_format(std::istream &,
		                                                      const file::pdb &,
		                                                      const file::pdb &,
		                                                      std::ostream &arg_stderr = std::cerr);

		alignment read_alignment_from_cath_ssap_legacy_format(std::istream &,
		                                                      const residue_name_vec &,
		                                                      const residue_name_vec &,
		                                                      std::ostream &arg_stderr = std::cerr);

		alignment read_alignment_from_cath_cora_legacy_format(std::istream &,
		                                                      const file::pdb_list &,
		                                                      std::ostream &arg_stderr = std::cerr);

		str_str_pair_vec read_ids_and_sequences_from_fasta(std::istream &);

		aln_posn_opt_vec align_sequence_to_amino_acids(const std::string &,
		                                               const amino_acid_vec &,
		                                               const std::string &,
		                                               std::ostream &arg_stderr = std::cerr);

		/// \todo Write code to make it easy to parse from proteins (and to optionally rescore whilst doing so)

		alignment read_alignment_from_fasta_file(const boost::filesystem::path &,
		                                         const file::pdb_list &,
		                                         std::ostream &arg_stderr = std::cerr);

		alignment read_alignment_from_fasta_file(const boost::filesystem::path &,
		                                         const file::pdb_list &,
		                                         const str_vec &,
		                                         std::ostream &arg_stderr = std::cerr);

		alignment read_alignment_from_fasta(std::istream &,
		                                    const file::pdb_list &,
		                                    std::ostream &arg_stderr = std::cerr);

		alignment read_alignment_from_fasta_file(const boost::filesystem::path &,
		                                         const protein_list &,
		                                         std::ostream &arg_stderr = std::cerr);

		alignment read_alignment_from_fasta_file(const boost::filesystem::path &,
		                                         const protein_list &,
		                                         const str_vec &,
		                                         std::ostream &arg_stderr = std::cerr);

		alignment read_alignment_from_fasta(std::istream &,
		                                    const protein_list &,
		                                    std::ostream &arg_stderr = std::cerr);

		alignment read_alignment_from_fasta(std::istream &,
		                                    const amino_acid_vec_vec &,
		                                    const str_vec &,
		                                    std::ostream &arg_stderr = std::cerr);

		aln_posn_opt search_for_residue_in_residue_names(const size_t &,
		                                                 const residue_name_vec &,
		                                                 const char &,
		                                                 const residue_name &,
		                                                 std::ostream &arg_stderr = std::cerr);

		void write_alignment_as_cath_ssap_legacy_format(const boost::filesystem::path &,
		                                                const alignment &,
		                                                const protein &,
		                                                const protein &);

		std::ostream & output_alignment_to_cath_ssap_legacy_format(std::ostream &,
		                                                           const alignment &,
		                                                           const protein &,
		                                                           const protein &);

		void write_alignment_as_fasta_alignment(const boost::filesystem::path &,
		                                        const alignment &,
		                                        const protein_list &);

		std::ostream & write_alignment_as_fasta_alignment(std::ostream &,
		                                                  const alignment &,
		                                                  const protein_list &);

		std::ostream & write_alignment_as_fasta_alignment(std::ostream &,
		                                                  const alignment &,
		                                                  const file::pdb_list &,
		                                                  const str_vec &);

		std::string alignment_as_fasta_string(const alignment &,
		                                      const protein_list &);

		std::string alignment_as_fasta_string(const alignment &,
		                                      const file::pdb_list &,
		                                      const str_vec &);
	} // namespace align
} // namespace cath

#endif
