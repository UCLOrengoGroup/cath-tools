/// \file
/// \brief The aligned_pair_score_value_list class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST_ALIGNED_PAIR_SCORE_VALUE_LIST_H
#define _CATH_TOOLS_SOURCE_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST_ALIGNED_PAIR_SCORE_VALUE_LIST_H

#include <boost/property_tree/ptree_fwd.hpp>

#include "score/aligned_pair_score/aligned_pair_score.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_list.hpp"

#include <vector>

namespace cath {
	namespace align {
		class alignment;
	} // namespace align

	namespace score {

		/// \brief TODOCUMENT
		class aligned_pair_score_value_list final {
		private:
			/// \brief TODOCUMENT
			score_value_vec         score_values;

			/// \brief TODOCUMENT
			aligned_pair_score_list scores;

		public:
			void add_score_and_value(const aligned_pair_score &,
			                         const score_value &);

			size_t size() const;
			score_value get_value_of_index(const size_t &) const;
			const aligned_pair_score & get_aligned_pair_score_of_index(const size_t &) const;
		};

		void warn_on_duplicate_human_friendly_names(const aligned_pair_score_value_list &);

		str_vec get_all_names(const aligned_pair_score_value_list &);

		void save_to_ptree(boost::property_tree::ptree &,
		                   const aligned_pair_score_value_list &);

		aligned_pair_score_value_list make_aligned_pair_score_value_list(const aligned_pair_score_list &,
		                                                                 const align::alignment &,
		                                                                 const protein &,
		                                                                 const protein &);

		std::ostream & operator<<(std::ostream &,
		                          const aligned_pair_score_value_list &);
	} // namespace score

} // namespace cath

#endif
