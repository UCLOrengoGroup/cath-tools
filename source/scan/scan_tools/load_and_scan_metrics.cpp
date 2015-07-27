/// \file
/// \brief The load_and_scan_metrics class definitions

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

#include "load_and_scan_metrics.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>
#include <boost/units/io.hpp>

#include "common/chrono/duration_to_seconds_string.h"
#include "common/file/open_fstream.h"
#include "common/type_aliases.h"

#include <fstream>

using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace cath::scan;
using namespace cath::scan::detail;
//using namespace std;

using boost::algorithm::join;
using boost::adaptors::transformed;
using boost::lexical_cast;
using boost::range::max_element;
using std::ofstream;
using std::string;

/// \brief TODOCUMENT
load_and_scan_metrics::load_and_scan_metrics(const hrc_duration &arg_load_files_durn, ///< TODOCUMENT
                                             const scan_metrics &arg_scan_metrics     ///< TODOCUMENT
                                             ) : load_files_durn  ( arg_load_files_durn ),
                                                 the_scan_metrics ( arg_scan_metrics    ) {
}

/// \brief TODOCUMENT
const hrc_duration & load_and_scan_metrics::get_load_files_durn() const {
	return load_files_durn;
}

/// \brief TODOCUMENT
const scan_metrics & load_and_scan_metrics::get_scan_metrics() const {
	return the_scan_metrics;
}

/// \brief TODOCUMENT
///
/// \relates load_and_scan_metrics
const durn_mem_pair & cath::scan::get_query_strucs_metrics(const load_and_scan_metrics &arg_load_and_scan_metrics ///< TODOCUMENT
                                                           ) {
	return arg_load_and_scan_metrics.get_scan_metrics().get_query_strucs_metrics();
}

/// \brief TODOCUMENT
///
/// \relates load_and_scan_metrics
const durn_mem_pair & cath::scan::get_query_index_metrics(const load_and_scan_metrics &arg_load_and_scan_metrics ///< TODOCUMENT
                                                          ) {
	return arg_load_and_scan_metrics.get_scan_metrics().get_query_index_metrics();
}

/// \brief TODOCUMENT
///
/// \relates load_and_scan_metrics
const durn_mem_pair & cath::scan::get_index_strucs_metrics(const load_and_scan_metrics &arg_load_and_scan_metrics ///< TODOCUMENT
                                                           ) {
	return arg_load_and_scan_metrics.get_scan_metrics().get_index_strucs_metrics();
}

/// \brief TODOCUMENT
///
/// \relates load_and_scan_metrics
const durn_mem_pair & cath::scan::get_index_index_metrics(const load_and_scan_metrics &arg_load_and_scan_metrics ///< TODOCUMENT
                                                          ) {
	return arg_load_and_scan_metrics.get_scan_metrics().get_index_index_metrics();
}

/// \brief TODOCUMENT
///
/// \relates load_and_scan_metrics
const hrc_duration & cath::scan::get_scan_durn(const load_and_scan_metrics &arg_load_and_scan_metrics ///< TODOCUMENT
                                               ) {
	return arg_load_and_scan_metrics.get_scan_metrics().get_scan_durn();
}

/// \brief TODOCUMENT
///
/// \relates load_and_scan_metrics
string cath::scan::to_markdown_string(const load_and_scan_metrics &arg_load_and_scan_metrics ///< TODOCUMENT
                                      ) {
	const auto query_strucs_metrics = get_query_strucs_metrics( arg_load_and_scan_metrics );
	const auto query_index_metrics  = get_query_index_metrics ( arg_load_and_scan_metrics );
	const auto index_strucs_metrics = get_index_strucs_metrics( arg_load_and_scan_metrics );
	const auto index_index_metrics  = get_index_index_metrics ( arg_load_and_scan_metrics );

	const auto &load_files_durn     = arg_load_and_scan_metrics.get_load_files_durn();
	const auto &query_strucs_durn   = query_strucs_metrics.first;
	const auto &query_index_durn    = query_index_metrics.first;
	const auto &index_strucs_durn   = index_strucs_metrics.first;
	const auto &index_index_durn    = index_index_metrics.first;
	const auto &scan_durn           = get_scan_durn           ( arg_load_and_scan_metrics );
	const auto total_durn           =   load_files_durn
	                                  + query_strucs_durn
	                                  + query_index_durn
	                                  + index_strucs_durn
	                                  + index_index_durn
	                                  + scan_durn;

	const auto &query_strucs_size   = query_strucs_metrics.second;
	const auto &query_index_size    = query_index_metrics.second;
	const auto &index_strucs_size   = index_strucs_metrics.second;
	const auto &index_index_size    = index_index_metrics.second;
	const auto  total_size          =   query_strucs_size
	                                  + query_index_size
	                                  + index_strucs_size
	                                  + index_index_size;

	const auto property_fields = str_str_str_str_tpl_vec{ {
		str_str_str_str_tpl{ "Task",                       "Duration",                                  "Rate",                                              "Memory Required"                         },
		str_str_str_str_tpl{ "Load files",                 durn_to_seconds_string( load_files_durn   ), durn_to_rate_per_second_string( load_files_durn   ), ""                                        },
		str_str_str_str_tpl{ "Build query structure data", durn_to_seconds_string( query_strucs_durn ), durn_to_rate_per_second_string( query_strucs_durn ), lexical_cast<string>( query_strucs_size ) },
		str_str_str_str_tpl{ "Build query index store",    durn_to_seconds_string( query_index_durn  ), durn_to_rate_per_second_string( query_index_durn  ), lexical_cast<string>( query_index_size  ) },
		str_str_str_str_tpl{ "Build match structure data", durn_to_seconds_string( index_strucs_durn ), durn_to_rate_per_second_string( index_strucs_durn ), lexical_cast<string>( index_strucs_size ) },
		str_str_str_str_tpl{ "Build match index store",    durn_to_seconds_string( index_index_durn  ), durn_to_rate_per_second_string( index_index_durn  ), lexical_cast<string>( index_index_size  ) },
		str_str_str_str_tpl{ "Build scan_duration",        durn_to_seconds_string( scan_durn         ), durn_to_rate_per_second_string( scan_durn         ), ""                                        },
		str_str_str_str_tpl{ "",                           "",                                          "",                                                  ""                                        },
		str_str_str_str_tpl{ "**Everything**",             durn_to_seconds_string( total_durn        ), durn_to_rate_per_second_string( total_durn        ), lexical_cast<string>( total_size        ) }
	} };


	return markdown_table( property_fields );
}

/// \brief TODOCUMENT
///
/// \relates load_and_scan_metrics
void cath::scan::to_markdown_file(const load_and_scan_metrics &arg_scan_metrics,    ///< TODOCUMENT
                                  const path                  &arg_markdown_outfile ///< TODOCUMENT
                                  ) {
	ofstream test_ofstream;
	open_ofstream( test_ofstream, arg_markdown_outfile );
	test_ofstream << to_markdown_string( arg_scan_metrics );
	test_ofstream.close();
}
