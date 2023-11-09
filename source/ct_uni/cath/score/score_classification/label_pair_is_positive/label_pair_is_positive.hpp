/// \file
/// \brief The label_pair_is_positive class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE_LABEL_PAIR_IS_POSITIVE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE_LABEL_PAIR_IS_POSITIVE_HPP

#include <filesystem>

#include "cath/common/type_aliases.hpp"

namespace cath::score {

	/// \brief Class to store pairs of labels along with whether they are positive (eg represent homologous domains) or negative
	class label_pair_is_positive final {
	private:
		/// \brief A map from the pairs of names to an is_positive bool
		str_str_pair_bool_map all_pairs;

	public:
		explicit label_pair_is_positive(str_str_pair_bool_map);

		[[nodiscard]] bool is_positive( const std::string &, const std::string & ) const;
	};

	label_pair_is_positive make_label_pair_is_positive(std::istream &);

	label_pair_is_positive make_label_pair_is_positive(const ::std::filesystem::path &);

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE_LABEL_PAIR_IS_POSITIVE_HPP
