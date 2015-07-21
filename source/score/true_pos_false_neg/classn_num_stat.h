/// \file
/// \brief The classn_num_stat class header

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

#ifndef CLASSN_NUM_STAT_H_INCLUDED
#define CLASSN_NUM_STAT_H_INCLUDED

#include "score/true_pos_false_neg/classn_outcome.h"
#include "score/true_pos_false_neg/classn_stat.h"
#include "score/true_pos_false_neg/true_false_pos_neg.h"

namespace cath {
	namespace score {

		/// \todo Concrete classn_stat for direct statistics of {TP, TN, FP, FN})
		///
		/// This allows these statistics to be simply defined as type aliases.
		template <classn_outcome N>
		class classn_num_stat final : public classn_stat {
			virtual size_rational do_calculate(const true_false_pos_neg &) const override final;

		public:
			virtual ~classn_num_stat() noexcept = default;
		};

		/// \brief Calculate the numerator and denominator of the rate and return a rational
		template <classn_outcome N>
		size_rational classn_num_stat<N>::do_calculate(const true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg from which to calculate the rate
		                                               ) const {
			const size_t &numerator = arg_true_false_pos_neg.get_num<N>();
			return size_rational( numerator, 1 );
		}

		using true_positive_stat  = classn_num_stat< classn_outcome::TRUE_POSITIVE  >;
		using true_negative_stat  = classn_num_stat< classn_outcome::TRUE_NEGATIVE  >;
		using false_positive_stat = classn_num_stat< classn_outcome::FALSE_POSITIVE >;
		using false_negative_stat = classn_num_stat< classn_outcome::FALSE_NEGATIVE >;
	}
}

#endif
