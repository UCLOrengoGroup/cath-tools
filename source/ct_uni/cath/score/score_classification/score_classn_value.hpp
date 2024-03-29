/// \file
/// \brief The score_classn_value class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_SCORE_CLASSN_VALUE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_SCORE_CLASSN_VALUE_HPP

#include <string>

namespace cath::score {

	/// \brief TODOCUMENT
	class score_classn_value final {
	private:
		/// \brief TODOCUMENT
		double score_value;

		/// \brief TODOCUMENT
		bool instance_is_positive;

		/// \brief TODOCUMENT
		std::string instance_label;

	public:
		score_classn_value(const double &,
		                   const bool &,
		                   std::string = {} );

		[[nodiscard]] const double &     get_score_value() const;
		[[nodiscard]] const bool &       get_instance_is_positive() const;
		[[nodiscard]] const std::string &get_instance_label() const;
	};

	std::ostream & operator<<(std::ostream &,
	                          const score_classn_value &);

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_SCORE_CLASSN_VALUE_HPP
