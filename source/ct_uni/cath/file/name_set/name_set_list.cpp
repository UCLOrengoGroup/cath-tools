/// \file
/// \brief The name_set_list class definitions

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

#include "name_set_list.hpp"

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/irange.hpp>

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::file;

using boost::adaptors::transformed;
using boost::algorithm::all_of;
using boost::algorithm::join;
using boost::none;
using boost::range::combine;
using std::min;
using std::ostream;
using std::string;

/// \brief Build a name_set_list from str_vec and a str_opt_vec
///
/// The resulting name_set_list will have the same size as prm_names_from_acq
///
/// \relates name_set_list
name_set_list cath::file::build_name_set_list(str_vec     prm_names_from_acq, ///< The names obtained from a pdbs_acquirer
                                              str_vec     prm_ids,            ///< Alternative IDs
                                              str_opt_vec prm_domains         ///< Regions for the strucs_context
                                              ) {
	// Convert ids to a new str_opt of the same size as prm_names_from_acq
	str_opt_vec ids;
	ids.reserve( prm_names_from_acq.size() );
	for (string &the_string : prm_ids) {
		ids.push_back( std::move( the_string ) );
	}
	ids.resize( prm_names_from_acq.size(), none );

	// Adjust prm_domains to the correct size (if necessary)
	prm_domains.resize( prm_names_from_acq.size(), none );

	// Build a name_set_vec from the data
	name_set_vec result;
	result.reserve( prm_names_from_acq.size() );
	for (boost::tuple<string &, str_opt &, str_opt &> &&x : combine( prm_names_from_acq, ids, prm_domains ) ) {
		result.emplace_back(
			std::move( x.get<0>() ),
			std::move( x.get<1>() ),
			std::move( x.get<2>() )
		);
	}

	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return name_set_list{ result };
}

/// \brief Generate a string describing the specified name_set_list
///
/// \relates name_set_list
string cath::file::to_string(const name_set_list &prm_name_set_list ///< The name_set_list to describe
                             ) {
	return "name_set_list[\n\t"
		+ join(
			prm_name_set_list
				| transformed( [] (const name_set &x) { return to_string( x ); } ),
			",\n\t"
		)
		+ "\n]";
}

/// \brief Insert a description of the specified name_set_list into the specified ostream
///
/// \relates name_set_list
ostream & cath::file::operator<<(ostream             &prm_os,           ///< The ostream into which the description should be inserted
                                 const name_set_list &prm_name_set_list ///< The name_set_list to describe
                                 ) {
	prm_os << to_string( prm_name_set_list );
	return prm_os;
}


/// \brief Return whether all the specified name_sets have specified_ids
///
/// \relates name_set_list
bool cath::file::all_have_specified_id(const name_set_list &prm_name_set_list ///< TODOCUMENT
                                       ) {
	return all_of(
		prm_name_set_list,
		[] (const name_set &x) { return static_cast<bool>( x.get_specified_id() ); }
	);
}

/// \brief Return whether all the specified name_sets have domain_name_from_regions
///
/// \relates name_set_list
bool cath::file::all_have_domain_name_from_regions(const name_set_list &prm_name_set_list ///< TODOCUMENT
                                                   ) {
	return all_of(
		prm_name_set_list,
		[] (const name_set &x) { return static_cast<bool>( x.get_domain_name_from_regions() ); }
	);
}

/// \brief Get a vector of the name_from_acq strings of the specified name_sets
///
/// \relates name_set_list
str_vec cath::file::get_names_from_acq(const name_set_list &prm_name_sets ///< The name_set_list to query
                                       ) {
	return transform_build<str_vec>(
		prm_name_sets,
		[] (const name_set &x) { return x.get_name_from_acq(); }
	);
}

/// \brief Get a vector of the name_from_acq strings of the specified name_sets
///
/// \relates name_set_list
str_vec cath::file::get_domain_or_specified_or_from_acq_names(const name_set_list &prm_name_sets ///< The name_set_list to query
                                                              ) {
	return transform_build<str_vec>(
		prm_name_sets,
		[] (const name_set &x) { return get_domain_or_specified_or_name_from_acq( x ); }
	);
}

/// \brief Get a vector of the names from the specified name_sets suitable for use in alignment HTML
///
/// \relates name_set_list
str_vec cath::file::get_alignment_html_names(const name_set_list &prm_name_sets ///< The name_set_list to query
                                             ) {
	return get_domain_or_specified_or_from_acq_names( prm_name_sets );
}

/// \brief Get a vector of the names from the specified name_sets suitable for use in generating multi_ssap_alignment file names
///
/// \relates name_set_list
str_vec cath::file::get_multi_ssap_alignment_file_names(const name_set_list &prm_name_sets ///< The name_set_list to query
                                                        ) {
	return get_domain_or_specified_or_from_acq_names( prm_name_sets );
}

/// \brief Get a vector of the names from the specified name_sets suitable for use in building protein_lists
///
/// \relates name_set_list
str_vec cath::file::get_protein_list_names(const name_set_list &prm_name_sets ///< The name_set_list to query
                                           ) {
	return get_domain_or_specified_or_from_acq_names( prm_name_sets );
}

/// \brief Get a vector of the names from the specified name_sets suitable for use in superposition JSON
///
/// \relates name_set_list
str_vec cath::file::get_supn_json_names(const name_set_list &prm_name_sets ///< The name_set_list to query
                                        ) {
	return get_names_from_acq( prm_name_sets );
}

/// \brief Get a vector of the names from the specified name_sets suitable for use in generating superposition pdb file names
///
/// \relates name_set_list
str_vec cath::file::get_supn_pdb_file_names(const name_set_list &prm_name_sets ///< The name_set_list to query
                                            ) {
	return get_domain_or_specified_or_from_acq_names( prm_name_sets );
}

/// \brief Get a vector of the names from the specified name_sets suitable for use in viewer (eg PyMOL) names
///
/// \relates name_set_list
str_vec cath::file::get_viewer_names(const name_set_list &prm_name_sets ///< The name_set_list to query
                                     ) {
	return get_domain_or_specified_or_from_acq_names( prm_name_sets );
}

/// \brief Add the specified specified IDs into the specified name_set_list
///
/// \relates name_set_list
void cath::file::add_specified_ids(name_set_list &prm_name_set_list, ///< The name_set_list to modify
                                   str_vec        prm_specified_ids  ///< The specified specified IDs to set
                                   ) {
	for (const size_t &index : indices( min( prm_name_set_list.size(), prm_specified_ids.size() ) ) ) {
		prm_name_set_list[ index ].set_specified_id(
			std::move( prm_specified_ids[ index ] )
		);
	}
}

/// \brief Copy the specified name_set_list, add the specified specified IDs and return
///
/// \relates name_set_list
name_set_list cath::file::add_specified_ids_copy(name_set_list prm_name_set_list, ///< The source name_set_list
                                                 str_vec       prm_specified_ids  ///< The specified specified IDs to set
                                                 ) {
	add_specified_ids( prm_name_set_list, std::move( prm_specified_ids ) );
	return prm_name_set_list;
}

/// \brief Add the specified domain names from regions into the specified name_set_list
///
/// \relates name_set_list
void cath::file::add_domain_names_from_regions(name_set_list &prm_name_set_list, ///< The name_set_list to modify
                                               str_opt_vec    prm_names          ///< The specified domain names to set
                                               ) {
	for (const size_t &index : indices( min( prm_name_set_list.size(), prm_names.size() ) ) ) {
		prm_name_set_list[ index ].set_domain_name_from_regions(
			std::move( prm_names[ index ] )
		);
	}
}

/// \brief Copy the specified name_set_list, add the specified domain names from regions and return
///
/// \relates name_set_list
name_set_list cath::file::add_domain_names_from_regions_copy(name_set_list prm_name_set_list, ///< The source name_set_list
                                                             str_opt_vec   prm_names          ///< The specified domain names to set
                                                             ) {
	add_domain_names_from_regions( prm_name_set_list, std::move( prm_names ) );
	return prm_name_set_list;
}
