/// \file
/// \brief The prc_scores_file class header

#ifndef PRC_SCORES_FILE_H_INCLUDED
#define PRC_SCORES_FILE_H_INCLUDED

#include <boost/filesystem/path.hpp>

#include "file/file_type_aliases.h"

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		class prc_scores_file final {
		private:
			prc_scores_file() = delete;

		public:
			static prc_scores_entry_vec remove_duplicates(const prc_scores_entry_vec &);

			static prc_scores_entry_vec parse_prc_scores_file(std::istream &);

			static prc_scores_entry_vec parse_prc_scores_file(const boost::filesystem::path &);

			static prc_scores_entry_vec parse_prc_scores_file_fancy(std::istream &);

			static prc_scores_entry_vec parse_prc_scores_file_fancy(const boost::filesystem::path &);
		};

	}
}

#endif
