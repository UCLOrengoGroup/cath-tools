/// \file
/// \brief The prc_scores_entry class header

#ifndef _CATH_TOOLS_SOURCE_FILE_PRC_SCORES_FILE_PRC_SCORES_ENTRY_H
#define _CATH_TOOLS_SOURCE_FILE_PRC_SCORES_FILE_PRC_SCORES_ENTRY_H

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
			prc_scores_entry(const std::string &,
			                 const size_t &,
			                 const size_t &,
			                 const size_t &,
			                 const size_t &,
			                 const std::string &,
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

#endif
