/// \file
/// \brief The prc_scores_entry class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PRC_SCORES_FILE_PRC_SCORES_ENTRY_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PRC_SCORES_FILE_PRC_SCORES_ENTRY_HPP

#include <boost/operators.hpp>

#include <string>

namespace cath {
	namespace file {

		/// \brief Represent data from a line of PRC data
		class prc_scores_entry final : private boost::equality_comparable<prc_scores_entry> {
		private:
			/// \brief Name of protein 1
			std::string name_1;

			/// \brief Start (of match?) on protein 1
			size_t start_1;

			/// \brief End (of match?) on protein 1
			size_t end_1;

			/// \brief Length of protein 1
			size_t length_1;

			/// \brief Number of this particular hit
			size_t hit_num;

			/// \brief Name of protein 2
			std::string name_2;

			/// \brief Start (of match?) on protein 2
			size_t start_2;

			/// \brief End (of match?) on protein 2
			size_t end_2;

			/// \brief Length of protein 2
			size_t length_2;

			/// \brief Simple score
			double simple;

			/// \brief Reverse score
			double reverse;

			/// \brief E-value
			double evalue;

		public:
			prc_scores_entry(std::string,
			                 const size_t &,
			                 const size_t &,
			                 const size_t &,
			                 const size_t &,
			                 std::string,
			                 const size_t &,
			                 const size_t &,
			                 const size_t &,
			                 const double &,
			                 const double &,
			                 const double &);

			const std::string & get_name_1() const;
			const size_t &      get_start_1() const;
			const size_t &      get_end_1() const;
			const size_t &      get_length_1() const;
			const size_t &      get_hit_num() const;
			const std::string & get_name_2() const;
			const size_t &      get_start_2() const;
			const size_t &      get_end_2() const;
			const size_t &      get_length_2() const;
			const double &      get_simple() const;
			const double &      get_reverse() const;
			const double &      get_evalue() const;
		};

		bool operator==(const prc_scores_entry &,
		                const prc_scores_entry &);

		prc_scores_entry prc_scores_entry_from_line(const std::string &);

		std::string to_string(const prc_scores_entry &);

		std::ostream & operator<<(std::ostream &,
		                          const prc_scores_entry &);

	} // namespace file
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PRC_SCORES_FILE_PRC_SCORES_ENTRY_HPP
