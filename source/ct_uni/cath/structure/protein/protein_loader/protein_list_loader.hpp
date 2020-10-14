/// \file
/// \brief The protein_list_loader class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_LOADER_PROTEIN_LIST_LOADER_HPP
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_LOADER_PROTEIN_LIST_LOADER_HPP

#include <boost/filesystem/path.hpp>

#include <iosfwd>
#include <utility>

#include "cath/common/type_aliases.hpp"
#include "cath/common/chrono/chrono_type_aliases.hpp"
#include "cath/common/clone/clone_ptr.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_source_file_set.hpp"

namespace cath { class protein_list; }

namespace cath {

	/// \brief Represent the details required to load a protein_list
	///	
	/// In the future, this can be changed to an ABC with concrete implementations
	/// for different approaches to loading proteins
	class protein_list_loader final {
	private:
		/// \brief The protein_source_file_set specifying which set of files should be used to build the protein
		const common::clone_ptr<protein_source_file_set> source_file_set_ptr;

		/// \brief The directory from which the files should be read
		const boost::filesystem::path data_dir;

		/// \brief The name of the proteins that are to be read from files
		const str_vec protein_names;

	public:
		protein_list_loader(const protein_source_file_set &,
		                    const boost::filesystem::path &,
		                    str_vec);

		std::pair<protein_list, hrc_duration> load_proteins(std::ostream &) const;
	};
} // namespace cath

#endif
