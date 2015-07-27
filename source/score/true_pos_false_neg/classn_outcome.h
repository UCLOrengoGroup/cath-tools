/// \file
/// \brief The classn_outcome header

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

#ifndef CLASSN_OUTCOME_H_INCLUDED
#define CLASSN_OUTCOME_H_INCLUDED

#include <iosfwd>

namespace cath {
	namespace score {

		/// \brief The outcome of trying to classify positives/negatives
		enum class classn_outcome {
			TRUE_POSITIVE,  ///< A correct    identification of a positive as a positive
			TRUE_NEGATIVE,  ///< A correct    identification of a negative as a negative
			FALSE_NEGATIVE, ///< An incorrect identification of a positive as a negative
			FALSE_POSITIVE  ///< An incorrect identification of a negative as a positive
		};

		classn_outcome outcome_of_correct_and_decision(const bool &,
		                                               const bool &);

		std::ostream & operator<<(std::ostream &,
		                          const classn_outcome &);

	}
}

#endif
