/// \file
/// \brief The classn_outcome definitions

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

#include "classn_outcome.hpp"

#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath::common;
using namespace ::cath::score;
using namespace ::std;

/// \brief Return the classn_outcome resulting from a bool correct answer and a bool decision
///
/// The result is "TRUE"     if the decision matches the answer, "FALSE" otherwise
/// The result is "POSITIVE" if the decision is true,            "FALSE" otherwise
classn_outcome cath::score::outcome_of_correct_and_decision(const bool &prm_correct_answer, ///< The correct answer for some boolean classification
                                                            const bool &prm_decision        ///< The decision of some boolean classification method
                                                            ) {
	return prm_correct_answer ? (prm_decision ? classn_outcome::TRUE_POSITIVE  : classn_outcome::FALSE_NEGATIVE )
	                          : (prm_decision ? classn_outcome::FALSE_POSITIVE : classn_outcome::TRUE_NEGATIVE  );
}

/// \brief TODOCUMENT
///
/// \relates classn_outcome
ostream & cath::score::operator<<(ostream              &prm_os,            ///< TODOCUMENT
                                  const classn_outcome &prm_classn_outcome ///< TODOCUMENT
                                  ) {
	switch ( prm_classn_outcome ) {
		case ( classn_outcome::TRUE_POSITIVE  ) : { prm_os << "TRUE_POSITIVE"  ; return prm_os ; }
		case ( classn_outcome::TRUE_NEGATIVE  ) : { prm_os << "TRUE_NEGATIVE"  ; return prm_os ; }
		case ( classn_outcome::FALSE_NEGATIVE ) : { prm_os << "FALSE_NEGATIVE" ; return prm_os ; }
		case ( classn_outcome::FALSE_POSITIVE ) : { prm_os << "FALSE_POSITIVE" ; return prm_os ; }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("protein_file_combn is not recognised"));
}
