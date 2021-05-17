/// \file
/// \brief The align_refining class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER_ALIGN_REFINING_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER_ALIGN_REFINING_HPP

#include <boost/any.hpp>

#include "cath/common/cpp20/make_array.hpp"
#include "cath/common/type_aliases.hpp"

#include <array>
#include <string>

namespace cath {
	namespace align {

		/// \brief Represent how much refining should be done to an alignment
		///        (typically as part of the alignment being acquired or immediately afterwards)
		enum class align_refining {
			NO,    ///< No refining should be performed
			LIGHT, ///< At most, light refining should be performed
			HEAVY  ///< Heavy, slow, expensive refining should be performed
		};

		/// \brief A constexpr list of all align_refinings
		[[maybe_unused]] constexpr auto all_align_refinings = common::make_array(
			align_refining::NO,
			align_refining::LIGHT,
			align_refining::HEAVY
		);

		std::string to_string(const align_refining &);

		std::ostream & operator<<(std::ostream &,
		                          const align_refining &);

		std::istream & operator>>(std::istream &,
		                          align_refining &);

		std::string description_of_align_refining(const align_refining &);

		void validate(boost::any &,
		              const str_vec &,
		              align_refining*,
		              int);

	} // namespace align
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER_ALIGN_REFINING_HPP
