/// \file
/// \brief The ssap_score_post_processing class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "ssap_score_post_processing.h"

#include <iostream>

using namespace cath::score;
using namespace std;

/// \brief TODOCUMENT
map<ssap_score_post_processing, string> name_of_ssap_score_post_processing::get() {
	return {
		{ ssap_score_post_processing::SIMPLE_NORMLS,          "ssap_score_post_processing::simple_normalise"           },
		{ ssap_score_post_processing::COMPLX_NORMLS,          "ssap_score_post_processing::complex_normalise"          },
		{ ssap_score_post_processing::SIMPLE_NORMLS_THEN_LOG, "ssap_score_post_processing::simple_normalise_then_log"  },
		{ ssap_score_post_processing::COMPLX_NORMLS_THEN_LOG, "ssap_score_post_processing::complex_normalise_then_log" },
		{ ssap_score_post_processing::LOG_THEN_SIMPLE_NORMLS, "ssap_score_post_processing::log_then_simple_normalise"  }
	};
}

/// \brief TODOCUMENT
///
/// \relates ssap_score_post_processing
ostream & cath::score::operator<<(ostream                          &arg_os,                        ///< TODOCUMENT
                                  const ssap_score_post_processing &arg_ssap_score_post_processing ///< TODOCUMENT
                                  ) {
	arg_os << name_of_ssap_score_post_processing::get().at( arg_ssap_score_post_processing );
	return arg_os;
}

/// \brief Whether the ssap_score_post_processing implies taking logs after normalisation
bool cath::score::has_post_log(const ssap_score_post_processing &arg_post_processing ///< The ssap_score_post_processing to assess
                               ) {
	return (
		arg_post_processing == ssap_score_post_processing::SIMPLE_NORMLS_THEN_LOG
		||
		arg_post_processing == ssap_score_post_processing::COMPLX_NORMLS_THEN_LOG
	);
}


/// \brief Whether the ssap_score_post_processing implies taking logs before normalisation
bool cath::score::has_pre_log(const ssap_score_post_processing &arg_post_processing ///< The ssap_score_post_processing to assess
                              ) {
	return (
		arg_post_processing == ssap_score_post_processing::LOG_THEN_SIMPLE_NORMLS
	);
}

/// \brief Whether the ssap_score_post_processing implies a simple normalisation
bool cath::score::normalisation_is_simple(const ssap_score_post_processing &arg_post_processing ///< The ssap_score_post_processing to assess
                                          ) {
	return (
		arg_post_processing == ssap_score_post_processing::SIMPLE_NORMLS
		||
		arg_post_processing == ssap_score_post_processing::SIMPLE_NORMLS_THEN_LOG
		||
		arg_post_processing == ssap_score_post_processing::LOG_THEN_SIMPLE_NORMLS
	);
}

