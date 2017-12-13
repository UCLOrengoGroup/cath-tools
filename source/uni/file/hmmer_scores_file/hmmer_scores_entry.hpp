/// \file
/// \brief The hmmer_scores_entry class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_FILE_HMMER_SCORES_FILE_HMMER_SCORES_ENTRY_HPP
#define _CATH_TOOLS_SOURCE_UNI_FILE_HMMER_SCORES_FILE_HMMER_SCORES_ENTRY_HPP

#include <boost/operators.hpp>

#include "file/hmmer_scores_file/hmmer_name_handling.hpp"

#include <string>

namespace cath {
	namespace file {
		namespace detail {
			std::string strip_header_name(const std::string &);
		} // namespace detail
	} // namespace file
} // namespace cath

namespace cath {
	namespace file {

		/// \brief Represent data from a line of HMMER data
		class hmmer_scores_entry final : private boost::equality_comparable<hmmer_scores_entry> {
		private:
			/// \brief Name of protein 1
			std::string name_1;

			/// \brief Accession of protein 1 (or usually "-")
			std::string accession_1;

			/// \brief Name of protein 2
			std::string name_2;

			/// \brief Accession of protein 2 (or usually "-")
			std::string accession_2;

			/// \brief E-value (full sequence)
			double full_sequence_evalue;

			/// \brief Score (full sequence):
			double full_sequence_score;

			/// \brief Bias (full sequence)
			double full_sequence_bias;

			/// \brief E-value (best 1 domain)
			double best_1_domain_evalue;

			/// \brief Score (best 1 domain)
			double best_1_domain_score;

			/// \brief Bias (best 1 domain)
			double best_1_domain_bias;

			/// \brief Expected number of domains
			double expected_num_doms;

			/// \brief "Number of discrete regions defined" (from HMMER documentation)
			size_t reg;

			/// \brief "Number of regions that appeared to be multidomain" (from HMMER documentation)
			size_t clu;

			/// \brief "For envelopes that were defined by stochastic traceback clustering, how many of them overlap
			///         other envelopes" (from HMMER documentation)
			size_t ov;

			/// \brief "The total number of envelopes" (from HMMER documentation)
			size_t env;

			/// \brief "Number of domains defined" (from HMMER documentation)
			size_t dom;

			/// \brief "Number of domains satisfying reporting threshold" (from HMMER documentation)
			size_t rep;

			/// \brief "Number of domains satisfying inclusion thresholds" (from HMMER documentation)
			size_t inc;

			/// \brief "The target’s description line, as free text" (from HMMER documentation)
			std::string description;

		public:
			hmmer_scores_entry(std::string,
			                   std::string,
			                   std::string,
			                   std::string,
			                   const double &,
			                   const double &,
			                   const double &,
			                   const double &,
			                   const double &,
			                   const double &,
			                   const double &,
			                   const size_t &,
			                   const size_t &,
			                   const size_t &,
			                   const size_t &,
			                   const size_t &,
			                   const size_t &,
			                   const size_t &,
			                   std::string);

			std::string get_name_1() const;
			std::string get_accession_1() const;
			std::string get_name_2() const;
			std::string get_accession_2() const;

			double      get_full_sequence_evalue() const;
			double      get_full_sequence_score() const;
			double      get_full_sequence_bias() const;
			double      get_best_1_domain_evalue() const;
			double      get_best_1_domain_score() const;
			double      get_best_1_domain_bias() const;
			double      get_expected_num_doms() const;
			size_t      get_reg() const;
			size_t      get_clu() const;
			size_t      get_ov() const;
			size_t      get_env() const;
			size_t      get_dom() const;
			size_t      get_rep() const;
			size_t      get_inc() const;

			std::string get_description() const;
		};

		bool operator==(const hmmer_scores_entry &,
		                const hmmer_scores_entry &);

		hmmer_scores_entry hmmer_scores_entry_from_line(const std::string &,
		                                                const hmmer_name_handling & = hmmer_name_handling::STRIP);

		std::string to_string(const hmmer_scores_entry &);

		std::ostream & operator<<(std::ostream &,
		                          const hmmer_scores_entry &);

	} // namespace file
} // namespace cath

#endif
