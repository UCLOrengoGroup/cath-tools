/// \file
/// \brief The score_value_list_json_outputter class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_OUTPUTTER_SCORE_VALUE_LIST_JSON_OUTPUTTER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_OUTPUTTER_SCORE_VALUE_LIST_JSON_OUTPUTTER_HPP

#include "cath/common/json_style.hpp"
#include "cath/common/type_aliases.hpp"

#include <iosfwd>

// clang-format off
namespace cath::score { class aligned_pair_score_value_list; }
// clang-format on

namespace cath::score {

	/// \brief Simple wrapper class for outputting an aligned_pair_score_value_list to an ostream in JSON format
	///
	/// Use like this:
	///
	///    cerr << score_value_list_json_outputter( my_score_value_list ) << endl;
	///
	/// This provides a convenient way for the user to choose a different format
	/// when outputting an aligned_pair_score_value_list to an ostream via the insertion operator
	class score_value_list_json_outputter final {
	private:
		/// \brief TODOCUMENT
		const aligned_pair_score_value_list &the_aligned_pair_score_value_list;

		/// \brief TODOCUMENT
		const common::json_style the_json_style = DEFAULT_JSON_STYLE;

	public:
		explicit score_value_list_json_outputter(const aligned_pair_score_value_list &,
		                                         const common::json_style & = DEFAULT_JSON_STYLE);

		/// \brief Default value for the style in which to output the JSON
		static constexpr common::json_style DEFAULT_JSON_STYLE = common::json_style::PRETTY;

		[[nodiscard]] const aligned_pair_score_value_list &get_aligned_pair_score_value_list() const;
		[[nodiscard]] const common::json_style &           get_json_style() const;
	};

	std::ostream & operator<<(std::ostream &,
	                          const score_value_list_json_outputter &);

} // namespace cath::score

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_OUTPUTTER_SCORE_VALUE_LIST_JSON_OUTPUTTER_HPP
