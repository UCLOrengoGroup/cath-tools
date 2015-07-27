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

#ifndef SCORE_VALUE_LIST_JSON_OUTPUTTER_H_INCLUDED
#define SCORE_VALUE_LIST_JSON_OUTPUTTER_H_INCLUDED

#include "common/type_aliases.h"

#include <iosfwd>

namespace cath {
	namespace score {
		class aligned_pair_score_value_list;
	}
}

namespace cath {
	namespace score {

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
			const bool pretty_print = true;

		public:
			score_value_list_json_outputter(const aligned_pair_score_value_list &,
			                                const bool &arg_pretty_print = true);

			const aligned_pair_score_value_list & get_aligned_pair_score_value_list() const;
			const bool & get_pretty_print() const;
		};

		std::ostream & operator<<(std::ostream &,
		                          const score_value_list_json_outputter &);

	}
}

#endif
