/// \file
/// \brief The cath_align_scorer class header

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

#ifndef _CATH_TOOLS_SOURCE_CATH_SCORE_ALIGN_CATH_ALIGN_SCORER_HPP
#define _CATH_TOOLS_SOURCE_CATH_SCORE_ALIGN_CATH_ALIGN_SCORER_HPP

#include "cath/common/type_aliases.hpp"

#include <iostream>

namespace cath { namespace opts { class cath_score_align_options; } }

namespace cath {

	/// \brief A class containing a public static method to do the real work of a call to the cath-score-align command
	///
	/// The details of the job to be done are passed in a cath_score_align_options object.
	class cath_align_scorer {
	public:
		cath_align_scorer() = delete;
		~cath_align_scorer() = delete;

		static void score(const opts::cath_score_align_options &,
		                  std::istream & = std::cin,
		                  std::ostream & = std::cout,
		                  std::ostream & = std::cerr);
	};

} // namespace cath

#endif
