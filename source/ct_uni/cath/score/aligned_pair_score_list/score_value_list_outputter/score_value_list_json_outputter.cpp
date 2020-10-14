/// \file
/// \brief The score_value_list_json_outputter class definitions

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

#include "score_value_list_json_outputter.hpp"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

#include "cath/score/aligned_pair_score_list/aligned_pair_score_value_list.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::score;
using namespace ::std;

using ::boost::property_tree::ptree;

constexpr json_style score_value_list_json_outputter::DEFAULT_JSON_STYLE;

/// \brief Ctor for score_value_list_json_outputter
score_value_list_json_outputter::score_value_list_json_outputter(const aligned_pair_score_value_list &prm_aligned_pair_score_value_list, ///< The alignment to be output
                                                                 const json_style                    &prm_json_style                     ///< The style in which the JSON should be written
                                                                 ) : the_aligned_pair_score_value_list ( prm_aligned_pair_score_value_list ),
                                                                     the_json_style                    ( prm_json_style                    ) {
}

/// \brief Getter for the const reference to the aligned_pair_score_value_list
const aligned_pair_score_value_list & score_value_list_json_outputter::get_aligned_pair_score_value_list() const {
	return the_aligned_pair_score_value_list;
}

/// \brief Getter for pretty_print
const json_style & score_value_list_json_outputter::get_json_style() const {
	return the_json_style;
}

/// \brief Output the aligned_pair_score_value_list to the ostream in JSON format
///
/// \relates score_value_list_json_outputter
ostream & cath::score::operator<<(ostream                               &prm_os,                             ///< The ostream to which the aligned_pair_score_value_list should be output
                                  const score_value_list_json_outputter &prm_score_value_list_json_outputter ///< A score_value_list_json_outputter that wraps the aligned_pair_score_value_list to be output
                                  ) {
	// Grab aligned_pair_score_value_list and then whether the JSON should be pretty printed
	const aligned_pair_score_value_list &the_aligned_pair_score_value_list = prm_score_value_list_json_outputter.get_aligned_pair_score_value_list();
	const json_style                    &the_json_style                    = prm_score_value_list_json_outputter.get_json_style();

	ptree temp_ptree;
	save_to_ptree( temp_ptree, the_aligned_pair_score_value_list );
	write_json( prm_os, temp_ptree, ( the_json_style == json_style::PRETTY ) );

	// Return the specified ostream
	return prm_os;
}
