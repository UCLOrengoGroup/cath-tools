/// \file
/// \brief The label_pair_is_positive class header

#ifndef _CATH_TOOLS_SOURCE_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE_LABEL_PAIR_IS_POSITIVE_H
#define _CATH_TOOLS_SOURCE_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE_LABEL_PAIR_IS_POSITIVE_H

#include <boost/filesystem/path.hpp>

#include "common/type_aliases.h"

namespace cath {
	namespace score {

		/// \brief Class to store pairs of labels along with whether they are positive (eg represent homologous domains) or negative
		class label_pair_is_positive final {
		private:
			/// \brief A map from the pairs of names to an is_positive bool
			str_str_pair_bool_map all_pairs;

		public:
			explicit label_pair_is_positive(const str_str_pair_bool_map &);

			bool is_positive(const std::string &,
			                 const std::string &) const;
		};

		label_pair_is_positive make_label_pair_is_positive(std::istream &);

		label_pair_is_positive make_label_pair_is_positive(const boost::filesystem::path &);

	} // namespace score
} // namespace cath

#endif
