/// \file
/// \brief The ssap_and_prc class header

#ifndef SSAP_AND_PRC_H_INCLUDED
#define SSAP_AND_PRC_H_INCLUDED

#include "file/prc_scores_file/prc_scores_entry.h"
#include "file/ssap_scores_file/ssap_scores_entry.h"

namespace cath { namespace file { class prc_scores_entry; } }
namespace cath { namespace file { class ssap_scores_entry; } }

namespace cath {
	namespace homcheck {

		/// \brief Represent a SSAP result and a PRC result for the same query/match pair
		class ssap_and_prc final {
		private:
			/// \brief The SSAP result
			file::ssap_scores_entry the_ssap;

			/// \brief The PRC result
			file::prc_scores_entry the_prc;

			/// \brief The magic function score calculated from the SSAP and PRC results
			///
			/// This gets populated by the ctor
			double magic_function_score;

		public:
			ssap_and_prc(const file::ssap_scores_entry &,
			             const file::prc_scores_entry &);

			const std::string & get_query_id() const;
			const std::string & get_match_id() const;

			const file::ssap_scores_entry & get_ssap() const;
			const file::prc_scores_entry & get_prc() const;
			const double & get_magic_function_score() const;
		};



		double magic_function(const file::ssap_scores_entry &,
		                      const file::prc_scores_entry &);

		const size_t & get_ssap_length_1(const ssap_and_prc &);
		const size_t & get_ssap_length_2(const ssap_and_prc &);
		const double & get_ssap_score(const ssap_and_prc &);
		const size_t & get_ssap_num_equivs(const ssap_and_prc &);
		const double & get_ssap_overlap_pc(const ssap_and_prc &);
		const double & get_ssap_seq_id_pc(const ssap_and_prc &);
		const double & get_ssap_rmsd(const ssap_and_prc &);
		const size_t & get_prc_start_1(const ssap_and_prc &);
		const size_t & get_prc_end_1(const ssap_and_prc &);
		const size_t & get_prc_length_1(const ssap_and_prc &);
		const size_t & get_prc_hit_num(const ssap_and_prc &);
		const size_t & get_prc_start_2(const ssap_and_prc &);
		const size_t & get_prc_end_2(const ssap_and_prc &);
		const size_t & get_prc_length_2(const ssap_and_prc &);
		const double & get_prc_simple(const ssap_and_prc &);
		const double & get_prc_reverse(const ssap_and_prc &);
		const double & get_prc_evalue(const ssap_and_prc &);

		std::string to_string(const ssap_and_prc &);
		std::ostream & operator<<(std::ostream &,
		                          const ssap_and_prc &);
	}
}

#endif
