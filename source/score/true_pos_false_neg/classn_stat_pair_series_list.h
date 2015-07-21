/// \file
/// \brief The classn_stat_pair_series_list class header

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

#ifndef CLASSN_STAT_PAIR_SERIES_LIST_H_INCLUDED
#define CLASSN_STAT_PAIR_SERIES_LIST_H_INCLUDED

#include "score/score_type_aliases.h"

namespace cath {
	namespace score {

		/// \brief TODOCUMENT
		class classn_stat_pair_series_list final {
		private:
			/// \brief TODOCUMENT
			classn_stat_pair_series_vec classn_stat_pair_serieses;

		public:
			/// \brief TODOCUMENT
			using const_iterator = classn_stat_pair_series_vec_citr;

			classn_stat_pair_series_list(const classn_stat_pair_series_vec &);

			bool empty() const;
			size_t size() const;

			const classn_stat_pair_series & operator[](const size_t &) const;
			const_iterator begin() const;
			const_iterator end() const;
		};

		const classn_stat_pair_series & classn_stat_pair_series_list_of_name(const classn_stat_pair_series_list &,
		                                                                     const std::string &);

	}
}

#endif
