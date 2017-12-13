/// \file
/// \brief The restrict_protein_source_file_set class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET_RESTRICT_PROTEIN_SOURCE_FILE_SET_HPP
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET_RESTRICT_PROTEIN_SOURCE_FILE_SET_HPP

#include "structure/protein/protein_source_file_set/protein_source_file_set.hpp"

namespace cath {

	/// \brief ABC for protein_source_file_set that provide their own do_read_and_restrict_files()
	///
	/// This disable the default do_read_and_restrict_files() and implements do_read_files() in terms of that
	class restrict_protein_source_file_set : public protein_source_file_set {
	private:
		protein do_read_files(const file::data_file_path_map &,
		                      const std::string &,
		                      std::ostream &) const final;

		/// \brief Virtual method with which each concrete restrict_protein_source_file_set must define how
		///        to do_read_files() and restrict the result to the specified regions.
		///
		/// This disables the default in protein_source_file_set
		virtual protein do_read_and_restrict_files(const file::data_file_path_map &,
		                                           const std::string &,
		                                           const chop::region_vec_opt &,
		                                           std::ostream &) const = 0;

	public:
		restrict_protein_source_file_set() = default;
		virtual ~restrict_protein_source_file_set() noexcept = default;
	};

} // namespace cath

#endif
