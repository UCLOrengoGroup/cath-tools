/// \file
/// \brief The ssap_scores class header

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

#ifndef SSAP_SCORES_H_INCLUDED
#define SSAP_SCORES_H_INCLUDED

#include <iosfwd>

namespace cath {

	/// \brief TODOCUMENT (Complex score measures)
	class ssap_scores final {
	private:
		/// \brief TODOCUMENT
		size_t num_aligned_pairs                    = 0;

		/// \brief TODOCUMENT
		double percentage_aligned_pairs_over_larger = 0.0;

		/// \brief TODOCUMENT
		double seq_id                               = 0.0;

		/// \brief TODOCUMENT
		double ssap_score_over_compared             = 0.0; ///< This doesn't appear to ever be used

		/// \brief TODOCUMENT
		double ssap_score_over_smaller              = 0.0;

		/// \brief TODOCUMENT
		double ssap_score_over_larger               = 0.0;

	public:
		void set_num_aligned_pairs(const size_t &);
		void set_percentage_aligned_pairs_over_larger(const double &);
		void set_seq_id(const double &);
		void set_ssap_score_over_compared(const double &);
		void set_ssap_score_over_smaller(const double &);
		void set_ssap_score_over_larger(const double &);

		size_t get_num_aligned_pairs() const;
		double get_percentage_aligned_pairs_over_larger() const;
		double get_seq_id() const;
		double get_ssap_score_over_compared() const;
		double get_ssap_score_over_smaller() const;
		double get_ssap_score_over_larger() const;
	};

	std::ostream & operator<<(std::ostream &,
	                          const ssap_scores &);

}
#endif
