/// \file
/// \brief The classn_stat class header

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

#ifndef _CATH_TOOLS_SOURCE_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_H
#define _CATH_TOOLS_SOURCE_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_H

#include <boost/rational.hpp>

#include "common/type_aliases.h"

namespace cath { using size_rational                             = boost::rational<size_t>;                              }
namespace cath { using size_rational_vec                         = std::vector<size_rational>;                           }
namespace cath { using size_rational_size_rational_pair          = std::pair<size_rational, size_rational>;              }
namespace cath { using size_rational_size_rational_pair_vec      = std::vector<size_rational_size_rational_pair>;        }
namespace cath { using size_rational_size_rational_pair_vec_itr  = size_rational_size_rational_pair_vec::iterator;       }
namespace cath { using size_rational_size_rational_pair_vec_citr = size_rational_size_rational_pair_vec::const_iterator; }

namespace cath { namespace common { class ratio; } }
namespace cath { namespace score { class true_false_pos_neg; } }

namespace cath {
	namespace score {

		/// \brief ABC for classification statistics (eg sensitivity, specificity etc) for true_false_pos_neg
		///
		/// \todo Since all of the statistics like sensitivity, specificity, precision etc
		///       are of the form A/(A+B) where A and B are one of {TP, TN, FP, FN}
		///       and since true_false_pos_neg provides a getter templated on classn_outcome,
		///       it'd be pretty easy to build a template that allows all of these statistics
		///       to be specified with two classn_outcome template values.
		///
		/// \todo Consider changing the return type to
		class classn_stat {
		private:
			/// \brief Pure virtual method with which each concrete classn_stat must define how
			/// to calculate the numerator and denominator of the statistic
			virtual size_rational do_calculate(const true_false_pos_neg &) const = 0;

			/// \brief TODOCUMENT
			virtual std::string do_get_name() const = 0;

		public:
			classn_stat() = default;
			virtual ~classn_stat() noexcept = default;

			classn_stat(const classn_stat &) = default;
			classn_stat(classn_stat &&) noexcept = default;
			classn_stat & operator=(const classn_stat &) = default;
			classn_stat & operator=(classn_stat &&) noexcept = default;

			size_rational calculate(const true_false_pos_neg &) const;
			std::string get_name() const;
		};

		double calculate_and_convert(const classn_stat &,
		                             const true_false_pos_neg &);

	} // namespace score
} // namespace cath

#endif
