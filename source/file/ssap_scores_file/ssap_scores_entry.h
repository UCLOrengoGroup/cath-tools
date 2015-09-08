/// \file
/// \brief The ssap_scores_entry class header

#ifndef SSAP_SCORES_ENTRY_H_INCLUDED
#define SSAP_SCORES_ENTRY_H_INCLUDED

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
			ssap_scores_entry(const std::string &,
			                  const std::string &,
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
	}
}

#endif
