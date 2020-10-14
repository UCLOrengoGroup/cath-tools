/// \file
/// \brief The residue_name_align_map class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL_RESIDUE_NAME_ALIGN_MAP_HPP
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL_RESIDUE_NAME_ALIGN_MAP_HPP

#include <iosfwd>
#include <map>
#include <vector>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace align {
		namespace detail {

			/// \brief Form and hold a map of the indices of each of the residue names in a list
			///
			/// This is a utility class to help implement alignment_coord_extractor.
			///
			/// It is basically a convenience wrapper for map<string, size_t> that:
			///  - auto-constructs from a vector of strings,
			///  - provides a single value_of_key method and
			///  - throws sensible exceptions in sensible places
			class residue_name_align_map final {
				/// \brief TODOCUMENT
				str_size_map index_of_residue_name;

			public:
				explicit residue_name_align_map(const str_vec &);

				bool contains_residue_name_string(const std::string &) const;
				size_t get_index_of_residue_name_string(const std::string &) const;

				str_vec get_residue_name_strings() const;
			};

			residue_name_align_map make_residue_name_align_map(const residue_name_vec &);
			bool contains_residue_name(const residue_name_align_map &,
			                           const residue_name &);
			size_t get_index_of_residue_name(const residue_name_align_map &,
			                                 const residue_name &);

			std::ostream & operator<<(std::ostream &,
			                          const residue_name_align_map &prm_residue_name_align_map);
		} // namespace detail
	} // namespace align
} // namespace cath

#endif
