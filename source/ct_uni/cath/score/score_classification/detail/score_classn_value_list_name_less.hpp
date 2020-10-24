/// \file
/// \brief The score_classn_value_list_name_less class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_DETAIL_SCORE_CLASSN_VALUE_LIST_NAME_LESS_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_DETAIL_SCORE_CLASSN_VALUE_LIST_NAME_LESS_HPP

#include <string>

namespace cath { namespace score { class score_classn_value_list; } }

namespace cath {
	namespace score {
		namespace detail {

			/// \brief TODOCUMENT
			class score_classn_value_list_name_less final {
			public:
				bool operator()(const score_classn_value_list &,
				                const score_classn_value_list &) const;

				bool operator()(const score_classn_value_list &,
				                const std::string &) const;
			};

		} // namespace detail
	} // namespace score
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_DETAIL_SCORE_CLASSN_VALUE_LIST_NAME_LESS_HPP
