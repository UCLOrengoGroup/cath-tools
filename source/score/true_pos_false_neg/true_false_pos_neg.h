/// \file
/// \brief The true_false_pos_neg class header

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

#ifndef TRUE_FALSE_POS_NEG_H_INCLUDED
#define TRUE_FALSE_POS_NEG_H_INCLUDED

#include <boost/operators.hpp>

#include "score/true_pos_false_neg/classn_outcome.h"

#include <cstddef>

namespace cath {
	namespace score {

		/// \brief Store the numbers of true/false positives/negatives associated with a classification attempt
		class true_false_pos_neg final : private boost::addable<true_false_pos_neg, classn_outcome> {
		private:
			/// \brief The number of true positives
			size_t num_true_positives  = 0;

			/// \brief The number of true negatives
			size_t num_true_negatives  = 0;

			/// \brief The number of false positives
			size_t num_false_positives = 0;

			/// \brief The number of false negatives
			size_t num_false_negatives = 0;

		public:
			true_false_pos_neg() = default;
			true_false_pos_neg(const size_t &,
			                   const size_t &,
			                   const size_t &,
			                   const size_t &);

			true_false_pos_neg & operator+=(const classn_outcome &);

			template <classn_outcome CLN>
			const size_t & get_num() const;

			template <classn_outcome CLN>
			void set_num(const size_t &);
		};

		/// \brief Getter for the number of some classn_outcome in a true_false_pos_neg
		template <classn_outcome CLN>
		inline const size_t & true_false_pos_neg::get_num() const {
			switch ( CLN ) {
				case ( classn_outcome::TRUE_POSITIVE  ) : { return num_true_positives;  }
				case ( classn_outcome::TRUE_NEGATIVE  ) : { return num_true_negatives;  }
				case ( classn_outcome::FALSE_POSITIVE ) : { return num_false_positives; }
				case ( classn_outcome::FALSE_NEGATIVE ) : { return num_false_negatives; }
			}
		}

		/// \brief Setter for the number of some classn_outcome in a true_false_pos_neg
		template <classn_outcome CLN>
		inline void true_false_pos_neg::set_num(const size_t &arg_number ///< The value to which to the number should be set
		                                        ) {
			switch ( CLN ) {
				case ( classn_outcome::TRUE_POSITIVE  ) : { num_true_positives  = arg_number; break; }
				case ( classn_outcome::TRUE_NEGATIVE  ) : { num_true_negatives  = arg_number; break; }
				case ( classn_outcome::FALSE_POSITIVE ) : { num_false_positives = arg_number; break; }
				case ( classn_outcome::FALSE_NEGATIVE ) : { num_false_negatives = arg_number; break; }
			}
		}

		/// \brief Incrementer for the number of some classn_outcome in the specified true_false_pos_neg
		///
		/// \relates true_false_pos_neg
		template <classn_outcome CLN>
		inline void increment_num(true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg to be updated
		                          ) {
			arg_true_false_pos_neg.set_num<CLN>(
				1 + arg_true_false_pos_neg.get_num<CLN>()
			);
		}

		const size_t & get_num_true_positives(const true_false_pos_neg &);
		const size_t & get_num_true_negatives(const true_false_pos_neg &);
		const size_t & get_num_false_positives(const true_false_pos_neg &);
		const size_t & get_num_false_negatives(const true_false_pos_neg &);

		void set_num_true_positives(true_false_pos_neg &, const size_t &);
		void set_num_true_negatives(true_false_pos_neg &, const size_t &);
		void set_num_false_positives(true_false_pos_neg &, const size_t &);
		void set_num_false_negatives(true_false_pos_neg &, const size_t &);

		void increment_num_true_positives(true_false_pos_neg &);
		void increment_num_true_negatives(true_false_pos_neg &);
		void increment_num_false_positives(true_false_pos_neg &);
		void increment_num_false_negatives(true_false_pos_neg &);

		void update_with_predicted_positives(true_false_pos_neg &,
		                                     const size_t &,
											 const size_t &);

		std::ostream & operator<<(std::ostream &,
		                          const true_false_pos_neg &);
	}
}

#endif
