/// \file
/// \brief The score_classn_value_better_value class header

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

#ifndef SCORE_CLASSN_VALUE_BETTER_VALUE_H_INCLUDED
#define SCORE_CLASSN_VALUE_BETTER_VALUE_H_INCLUDED

namespace cath { namespace score { class score_classn_value; } }

namespace cath {
	namespace score {

		/// \brief Functor to return whether the first score_classn_value is "better" than the second
		///
		/// Better is defined based on the score_values using the direction specified on construction.
		/// If higher_is_better is false, this acts like a standard less-than functor on the score_values
		/// of the score_classn_values.
		///
		/// Used in a standard sort function, this sorts the best entries to the front.
		class score_classn_value_better_value final {
		private:
			/// \brief Whether or not a higher score_value is "better" and should be treated as less-than
			bool higher_is_better;

		public:
			/// \brief TODOCUMENT
			using first_argument_type  = const score_classn_value &;

			/// \brief TODOCUMENT
			using second_argument_type = const score_classn_value &;

			/// \brief TODOCUMENT
			using result_type          = bool;

			score_classn_value_better_value(const bool &);

			const bool & get_higher_is_better() const;

			bool operator()(const score_classn_value &,
			                const score_classn_value &);
		};

		double get_worst_possible_score(const score_classn_value_better_value &);
		double get_worst_possible_score(const bool &);

	}
}

#endif
