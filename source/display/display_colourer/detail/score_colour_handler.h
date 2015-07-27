/// \file
/// \brief The score_colour_handler class header

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

#ifndef SCORE_COLOUR_HANDLER_H_INCLUDED
#define SCORE_COLOUR_HANDLER_H_INCLUDED

#include "common/type_aliases.h"

//#include "display/display_colour/display_colour_gradient.h"
//#include "display/display_colourer/display_colourer.h"

namespace cath {
	namespace align {
		class alignment;
	}
	class display_colour;
}

namespace cath {
	namespace detail {

		/// \brief TODOCUMENT
		class score_colour_handler final {
		private:
			/// \brief TODOCUMENT
			bool show_scores_if_present;

			/// \brief TODOCUMENT
			bool scores_to_equivs;

			/// \brief TODOCUMENT
			bool normalise_scores;

		public:
			score_colour_handler(const bool &,
			                     const bool &,
			                     const bool &);

			bool get_show_scores_if_present() const;
			bool get_scores_to_equivs() const;
			bool get_normalise_scores() const;

			float_score_type get_score_of_postion(const align::alignment &,
			                                      const size_t &,
			                                      const size_t &) const;
		};

		void score_colour(const score_colour_handler &,
		                  const align::alignment &,
		                  const size_t &,
		                  const size_t &,
		                  display_colour &);

		display_colour score_colour_copy(const score_colour_handler &,
		                                 const align::alignment &,
		                                 const size_t &,
		                                 const size_t &,
		                                 const display_colour);
	}
}

#endif
