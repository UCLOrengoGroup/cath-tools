/// \file
/// \brief The hmmer_scores_entry class header

#ifndef HMMER_SCORES_ENTRY_H_INCLUDED
#define HMMER_SCORES_ENTRY_H_INCLUDED

#include <boost/operators.hpp>

#include "file/hmmer_scores_file/hmmer_name_handling.h"

#include <string>

namespace cath {
	namespace file {
		namespace detail {
			std::string strip_header_name(const std::string &);
		}
	}
}

namespace cath {
	namespace file {

		/// \brief Represent data from a line of HMMer data
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

			/// \brief "Number of discrete regions defined" (from HMMer documentation)
			size_t reg;

			/// \brief "Number of regions that appeared to be multidomain" (from HMMer documentation)
			size_t clu;

			/// \brief "For envelopes that were defined by stochastic traceback clustering, how many of them overlap
			///         other envelopes" (from HMMer documentation)
			size_t ov;

			/// \brief "The total number of envelopes" (from HMMer documentation)
			size_t env;

			/// \brief "Number of domains defined" (from HMMer documentation)
			size_t dom;

			/// \brief "Number of domains satisfying reporting threshold" (from HMMer documentation)
			size_t rep;

			/// \brief "Number of domains satisfying inclusion thresholds" (from HMMer documentation)
			size_t inc;

			/// \brief "The target’s description line, as free text" (from HMMer documentation)
			std::string description;

		public:
			hmmer_scores_entry(const std::string &,
			                   const std::string &,
			                   const std::string &,
			                   const std::string &,
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
			                   const std::string &);

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
		                                                const hmmer_name_handling &arg_hmmer_name_handling = hmmer_name_handling::STRIP);

		std::string to_string(const hmmer_scores_entry &);

		std::ostream & operator<<(std::ostream &,
		                          const hmmer_scores_entry &);

	}
}

#endif
