/// \file
/// \brief The hmmer_scores_file class header

#ifndef HMMER_SCORES_FILE_H_INCLUDED
#define HMMER_SCORES_FILE_H_INCLUDED

#include <boost/filesystem/path.hpp>

#include "file/file_type_aliases.h"
#include "file/hmmer_scores_file/hmmer_name_handling.h"

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		class hmmer_scores_file final {
		private:
			hmmer_scores_file() = delete;

		public:
			static hmmer_scores_entry_vec remove_duplicates(const hmmer_scores_entry_vec &);

			static hmmer_scores_entry_vec parse_hmmer_scores_file(std::istream &,
			                                                      const hmmer_name_handling &arg_hmmer_name_handling = hmmer_name_handling::STRIP);

			static hmmer_scores_entry_vec parse_hmmer_scores_file(const boost::filesystem::path &,
			                                                      const hmmer_name_handling &arg_hmmer_name_handling = hmmer_name_handling::STRIP);

//			static hmmer_scores_entry_vec parse_hmmer_scores_file_fancy(std::istream &);

//			static hmmer_scores_entry_vec parse_hmmer_scores_file_fancy(const boost::filesystem::path &);
		};

	}
}

#endif
