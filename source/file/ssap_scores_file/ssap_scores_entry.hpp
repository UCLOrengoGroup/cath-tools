/// \file
/// \brief The ssap_scores_entry class header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_SSAP_SCORES_FILE_SSAP_SCORES_ENTRY_H
#define _CATH_TOOLS_SOURCE_FILE_SSAP_SCORES_FILE_SSAP_SCORES_ENTRY_H

#include <boost/operators.hpp>

#include <string>

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		///
		/// \todo Eradicate the substantial redundancy between this and cath::ssap_scores
		class ssap_scores_entry final : private boost::equality_comparable<ssap_scores_entry> {
		private:
			/// \brief Name of protein 1
			std::string name_1;

			/// \brief Name of protein 2
			std::string name_2;

			/// \brief Length of protein 1
			size_t      length_1;

			/// \brief Length of protein 2
			size_t      length_2;

			/// \brief SSAP score for structural comparison
			double      ssap_score;

			/// \brief Number of equivalent/aligned residues
			size_t      num_equivs;

			/// \brief Percentage overlap  (100% x overlap /length of largest)
			double      overlap_pc;

			/// \brief Percentage identity (100% x identity/length of smallest)
			double      seq_id_pc;

			/// \brief RMSD of superposed structures
			double      rmsd;

		public:
			ssap_scores_entry(std::string,
			                  std::string,
			                  const size_t &,
			                  const size_t &,
			                  const double &,
			                  const size_t &,
			                  const double &,
			                  const double &,
			                  const double &);

			const std::string & get_name_1() const;
			const std::string & get_name_2() const;
			const size_t & get_length_1() const;
			const size_t & get_length_2() const;
			const double & get_ssap_score() const;
			const size_t & get_num_equivs() const;
			const double & get_overlap_pc() const;
			const double & get_seq_id_pc() const;
			const double & get_rmsd() const;
		};

		bool operator==(const ssap_scores_entry &,
		                const ssap_scores_entry &);

		ssap_scores_entry ssap_scores_entry_from_line(const std::string &);

		std::string to_string(const ssap_scores_entry &);

		std::ostream & operator<<(std::ostream &,
		                          const ssap_scores_entry &);
	} // namespace file
} // namespace cath

#endif
