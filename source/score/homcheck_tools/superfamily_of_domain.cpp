/// \file
/// \brief The superfamily_of_domain class definitions

#include "superfamily_of_domain.h"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/file/open_fstream.h"
#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"

#include <fstream>
#include <regex>

using namespace cath::common;
using namespace cath::homcheck;
using namespace std;

using boost::algorithm::any_of;
using boost::algorithm::is_space;
using boost::filesystem::path;
using boost::token_compress_on;

/// \brief The string to use in between the fold of a new superfamily and the ID of the domain for which it's being created
const string superfamily_of_domain::NEW_SF_CORE_STRING = ".new_sf_in_fold_of_";

/// \brief Simple predicate to return whether a string is a valid superfamily_id
bool cath::homcheck::detail::is_valid_superfamily_id(const string &arg_superfamily_id_string ///< The string to check
                                                     ) {
    return regex_search( arg_superfamily_id_string, regex{ R"(^\d+\.\d+\.\d+\.\d+$)" } );
}

/// \brief Extract the fold from the specified superfamily ID string
///
/// It's assumed that the input is a valid superfamily ID or is constructed by a previous call to fold_of_superfamily_id()
string cath::homcheck::detail::fold_of_superfamily_id(const string &arg_superfamily_id ///< The superfamily ID from which the fold should be extracted
                                                      ) {
	return regex_replace( arg_superfamily_id, regex{ R"(\.[^\.]+$)" }, "" );
}

/// \brief Ctor from a vector<pair<string, string>> where each pair contains domain ID and the corresponding superfamily ID
superfamily_of_domain::superfamily_of_domain(const str_str_pair_vec &sf_of_dom ///< The domain ID -> superfamily ID data from which this superfamily_of_domain should be constructed
                                             ) : sf_of_dom( std::begin( sf_of_dom ), std::cend( sf_of_dom ) ) {
	for (const auto &x: sf_of_dom) {
		if (! detail::is_valid_superfamily_id( x.second ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to construct superfamily_of_domain with invalid superfamily ID "+ x.second));
		}
	}
}

/// \brief Getter for the size (ie number of domains for which superfamily information is stored)
size_t superfamily_of_domain::size() const {
	return sf_of_dom.size();
}

/// \brief Return whether or not this superfamily_of_domain has superfamily information for the specified domain ID
bool superfamily_of_domain::has_superfamily_of_domain(const string &arg_domain_id ///< The domain ID to query
                                                      ) const {
	return sf_of_dom.count( arg_domain_id ) > 0;
}

/// \brief Return whether or not this superfamily_of_domain has superfamily information for the specified domain ID
///
/// \pre `this->has_superfamily_of_domain( arg_domain_id )`, else and invalid_argument_exception will be thrown
const string & superfamily_of_domain::get_superfamily_of_domain(const string &arg_domain_id ///< The domain ID to query
                                                                ) const {
	const auto find_itr = sf_of_dom.find( arg_domain_id );
	if ( find_itr == std::cend( sf_of_dom ) ) {
		BOOST_THROW_EXCEPTION(
			invalid_argument_exception("Unable to find any entry in superfamily_of_domain for domain ID \""
			+ arg_domain_id
			+ "\""
		));
	}
	return find_itr->second;
}

/// \brief Add a new entry representing a new superfamily for the specified domain in the same fold as the second specified domain
///
/// \pre ` ! this->has_superfamily_of_domain( arg_new_domain_id   )`, else and invalid_argument_exception will be thrown
///          this->has_superfamily_of_domain( arg_match_domain_id )`, else and invalid_argument_exception will be thrown
void superfamily_of_domain::add_domain_in_new_sf_in_fold_of_domain(const string &arg_new_domain_id,  ///< The new domain for which an entry should be created
                                                                   const string &arg_match_domain_id ///< The existing domain in whose fold the new domain's new superfamily should be created
                                                                   ) {
	const string &superfamily_id = get_superfamily_of_domain( arg_match_domain_id );
	const string  fold_id        = detail::fold_of_superfamily_id( superfamily_id );

	if ( has_superfamily_of_domain( arg_new_domain_id) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot add entry for domain with ID \""
			+ arg_new_domain_id
			+ "\" into superfamily_of_domain because there it already has an entry"
		));
	}
	sf_of_dom.emplace( arg_new_domain_id, fold_id + NEW_SF_CORE_STRING + arg_match_domain_id );
}

/// \brief Parse the superfamily_of_domain information from the specified istream
///
/// \relates superfamily_of_domain
superfamily_of_domain cath::homcheck::parse_superfamily_of_domain(istream &arg_sf_of_dom_istream ///< The istream from which the superfamily_of_domain information should be parsed
                                                                  ) {
	str_str_pair_vec line_string_pairs;
	string line_string;
	while ( getline( arg_sf_of_dom_istream, line_string ) ) {
		const auto line_parts = split_build<str_vec>( line_string, is_space(), token_compress_on );
		if ( ! line_parts.empty() ) {
			if ( line_parts.size() != 2 ) {
				BOOST_THROW_EXCEPTION(runtime_error_exception(
					"Cannot parse superfamily_of_domain data because line \""
					+ line_string
					+ "\" does not contain two parts."
				));
			}
			const string &domain_id      = line_parts[ 0 ];
			const string &superfamily_id = line_parts[ 1 ];

			if ( ! detail::is_valid_superfamily_id( superfamily_id ) ) {
				BOOST_THROW_EXCEPTION(runtime_error_exception(
					"Second entry in parsed line from superfamily_of_domain data is \""
					+ superfamily_id
					+ "\", which isn't a valid superfamily ID"
				));
			}
			if ( any_of( line_string_pairs, [&] (const str_str_pair &x) { return ( x.first == domain_id ); } ) ) {
				BOOST_THROW_EXCEPTION(runtime_error_exception(
					"Duplicated domain_id \""
					+ domain_id
					+ "\" in superfamily_of_domain data"
				));
			}

			line_string_pairs.emplace_back( domain_id, superfamily_id );
		}
	}
	return { line_string_pairs };
}

/// \brief Parse the superfamily_of_domain information from the specified file
///
/// \relates superfamily_of_domain
superfamily_of_domain cath::homcheck::parse_superfamily_of_domain(const path &arg_sf_of_dom_file ///< The file from which the superfamily_of_domain information should be parsed
                                                                  ) {
	ifstream input_ifstream;
	open_ifstream( input_ifstream, arg_sf_of_dom_file );
	const auto sf_of_dom = parse_superfamily_of_domain( arg_sf_of_dom_file );
	input_ifstream.close();
	return sf_of_dom;
}

/// \brief Parse the superfamily_of_domain information from the specified string
///
/// \relates superfamily_of_domain
superfamily_of_domain cath::homcheck::parse_superfamily_of_domain(const string &arg_sf_of_dom_text ///< The string from which the superfamily_of_domain information should be parsed
                                                                  ) {
	istringstream input_ss{  arg_sf_of_dom_text};
	return parse_superfamily_of_domain( input_ss );
}
