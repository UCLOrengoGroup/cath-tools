/// \file
/// \brief The classn_stat class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "classn_stat.h"

#include "exception/invalid_argument_exception.h"

using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::rational_cast;

/// \brief NVI wrapper to the pure-virtual do_calculate with some post-processing
size_rational classn_stat::calculate(const true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg from which the statistic should be calculated
                                     ) const {
	// Call the pure-virtual do_calculate to get rational and
	// check the denominator isn't zero (else throw) before returning it
	const size_rational ratio = do_calculate( arg_true_false_pos_neg );
	if ( ratio.denominator() == 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot calculate classn_stat of true_false_pos_neg for which the denominator is zero"));
	}
	return ratio;
}

/// \brief TODOCUMENT
string classn_stat::get_name() const {
	return do_get_name();
}

/// \brief TODOCUMENT
double cath::score::calculate_and_convert(const classn_stat        &arg_classn_stat, ///< TODOCUMENT
                                          const true_false_pos_neg &arg_tfpn         ///< TODOCUMENT
                                          ) {
	return rational_cast<double>( arg_classn_stat.calculate( arg_tfpn ) );
}
