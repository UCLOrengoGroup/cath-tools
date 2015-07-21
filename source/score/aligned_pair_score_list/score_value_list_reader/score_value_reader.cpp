/// \file
/// \brief The score_value_reader class definitions

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

#include <boost/filesystem/path.hpp>
#include <boost/log/trivial.hpp>
//#include <boost/lexical_cast.hpp>
#include <boost/logic/tribool_io.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "common/file/open_fstream.h"
#include "score/aligned_pair_score/pseudo_string_score.h"
#include "score/aligned_pair_score_list/aligned_pair_score_value_list.h"
#include "score/aligned_pair_score_list/score_value_list_reader/score_value_reader.h"

using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::filesystem::path;
//using boost::lexical_cast;
using boost::logic::tribool;
using boost::property_tree::json_parser::read_json;
using boost::property_tree::ptree;

/// \brief TODOCUMENT
aligned_pair_score_value_list score_value_reader::read_aligned_pair_score_list_from_property_tree(const ptree &arg_property_tree ///< TODOCUMENT
                                                                                                  ) {
	aligned_pair_score_value_list results;
	for (const auto &the_score : arg_property_tree.get_child( "scores" ) ) {
		const auto score_type           = the_score.second.get<string>( "score_type"       );
		const auto score_value          = the_score.second.get<double>( "score_value"      );
		const auto higher_is_better_str = the_score.second.get<string>( "higher_is_better" );

		istringstream higher_is_better_ss( higher_is_better_str );
		tribool higher_is_better;
		higher_is_better_ss >> boolalpha >> higher_is_better;

		results.add_score_and_value(
			pseudo_string_score( score_type, higher_is_better ),
			score_value
		);
	}
	warn_on_duplicate_human_friendly_names( results );
	return results;
}

/// \brief TODOCUMENT
aligned_pair_score_value_list score_value_reader::read(istream &arg_istream ///< TODOCUMENT
                                                       ) {
	ptree temp_ptree;
	read_json( arg_istream, temp_ptree );
	return read_aligned_pair_score_list_from_property_tree( temp_ptree );
}

/// \brief TODOCUMENT
aligned_pair_score_value_list score_value_reader::read(const path &arg_input_file ///< TODOCUMENT
                                                       ) {
	BOOST_LOG_TRIVIAL( debug ) << "Reading aligned_pair_score_value_list from file " << arg_input_file;
	ifstream score_value_ifstream;
	open_ifstream( score_value_ifstream, arg_input_file );
	const aligned_pair_score_value_list result = read( score_value_ifstream );
	score_value_ifstream.close();
	return result;
}
