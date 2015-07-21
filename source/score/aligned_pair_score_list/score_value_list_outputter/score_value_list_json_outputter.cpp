/// \file
/// \brief The score_value_list_json_outputter class definitions

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

#include "score_value_list_json_outputter.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

#include "score/aligned_pair_score_list/aligned_pair_score_value_list.h"

//using namespace boost;
using namespace boost::property_tree;
using namespace cath;
using namespace cath::score;
using namespace std;

/// \brief Ctor for score_value_list_json_outputter
score_value_list_json_outputter::score_value_list_json_outputter(const aligned_pair_score_value_list &arg_aligned_pair_score_value_list, ///< The alignment to be output
                                                                 const bool                          &arg_pretty_print                   ///< TODOCUMENT
                                                                 ) : the_aligned_pair_score_value_list ( arg_aligned_pair_score_value_list ),
                                                                     pretty_print                      ( arg_pretty_print                  ) {
}

/// \brief Getter for the const reference to the aligned_pair_score_value_list
const aligned_pair_score_value_list & score_value_list_json_outputter::get_aligned_pair_score_value_list() const {
	return the_aligned_pair_score_value_list;
}

/// \brief Getter for pretty_print
const bool & score_value_list_json_outputter::get_pretty_print() const {
	return pretty_print;
}

/// \brief Output the aligned_pair_score_value_list to the ostream in JSON format
///
/// \relates score_value_list_json_outputter
ostream & cath::score::operator<<(ostream                               &arg_os,                             ///< The ostream to which the aligned_pair_score_value_list should be output
                                  const score_value_list_json_outputter &arg_score_value_list_json_outputter ///< A score_value_list_json_outputter that wraps the aligned_pair_score_value_list to be output
                                  ) {
	// Grab aligned_pair_score_value_list and then whether the JSON should be pretty printed
	const aligned_pair_score_value_list &the_aligned_pair_score_value_list = arg_score_value_list_json_outputter.get_aligned_pair_score_value_list();
	const bool                          &pretty_print                      = arg_score_value_list_json_outputter.get_pretty_print();

	ptree temp_ptree;
	save_to_ptree( temp_ptree, the_aligned_pair_score_value_list );
	write_json( arg_os, temp_ptree, pretty_print );

	// Return the specified ostream
	return arg_os;
}
