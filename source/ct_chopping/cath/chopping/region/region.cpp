/// \file
/// \brief The region class definitions

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

#include "region.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/program_options.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "cath/biocore/residue_id.hpp"
#include "cath/chopping/chopping_format/sillitoe_chopping_format.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/optional/make_optional_if.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::common;

using ::boost::adaptors::transformed;
using ::boost::algorithm::join;
using ::boost::any;
using ::boost::program_options::invalid_option_value;
using ::boost::program_options::validators::get_single_string;
using ::std::make_optional;
using ::std::make_pair;
using ::std::ostream;
using ::std::string;

/// \brief TODOCUMENT
void region::sanity_check() const {
	const residue_locating start_residue_locating = get_residue_locating( get_start_residue() );
	const residue_locating stop_residue_locating  = get_residue_locating( get_stop_residue()  );
	if ( start_residue_locating != stop_residue_locating ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to construct region with inconsistent approaches to identifying residues (names and/or indices)"));
	}

	if ( has_indices( *this ) ) {
		if ( get_start_index( *this ) > get_stop_index( *this ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct a region with start residue after the stop residue"));
		}
	}
}

/// \brief Ctor for a whole-chain region
region::region(const chain_label &prm_chain_label ///< The chain_label for the whole-chain region to construct
               ) : the_chain_label( prm_chain_label ) {
}

/// \brief Ctor for region
region::region(const chain_label  &prm_chain_label,        ///<TODOCUMENT
               const residue_name &prm_start_residue_name, ///<TODOCUMENT
               const residue_name &prm_stop_residue_name   ///<TODOCUMENT
               ) : the_chain_label( prm_chain_label        ),
                   residues       { make_pair(
                   	residue_location{ prm_start_residue_name },
                   	residue_location{ prm_stop_residue_name  }
                   ) } {
	sanity_check();
}

/// \brief Ctor for region
region::region(const chain_label  &prm_chain_label,        ///<TODOCUMENT
               const residue_name &prm_start_residue_name, ///<TODOCUMENT
               const size_t       &prm_start_index,        ///<TODOCUMENT
               const residue_name &prm_stop_residue_name,  ///<TODOCUMENT
               const size_t       &prm_stop_index          ///<TODOCUMENT
               ) : the_chain_label( prm_chain_label                         ),
                   residues       { make_pair(
                   	residue_location{ prm_start_residue_name, prm_start_index },
                   	residue_location{ prm_stop_residue_name,  prm_stop_index  }
                   ) } {
	sanity_check();
}

/// \brief Ctor for region
region::region(const size_t &prm_start_index, ///<TODOCUMENT
               const size_t &prm_stop_index   ///<TODOCUMENT
               ) : residues{ make_pair(
                   	residue_location{ prm_start_index },
                   	residue_location{ prm_stop_index  }
                   ) } {
	sanity_check();
}

/// \brief TODOCUMENT
const chain_label_opt & region::get_opt_chain_label() const {
	return the_chain_label;
}

/// \brief Return whether this region specified the start and stop (or is a whole-chain region)
bool region::has_starts_stops() const {
	return static_cast<bool>( residues );
}

/// \brief TODOCUMENT
const residue_location & region::get_start_residue() const {
	return residues->first;
}

/// \brief TODOCUMENT
const residue_location & region::get_stop_residue() const {
	return residues->second;
}

/// \brief Non-member equality operator for region
///
/// \relates region
bool cath::chop::operator==(const region &prm_lhs, ///< The first  region to compare
                            const region &prm_rhs  ///< The second region to compare
                            ) {
	return (
		( prm_lhs.get_opt_chain_label() == prm_rhs.get_opt_chain_label() )
		&&
		( prm_lhs.has_starts_stops()    == prm_rhs.has_starts_stops() )
		&&
		(
			prm_lhs.has_starts_stops()
				? (
					prm_lhs.get_start_residue()   == prm_rhs.get_start_residue()
					&&
					prm_lhs.get_stop_residue()    == prm_rhs.get_stop_residue()
				)
				: true
		)
	);
}

/// \brief TODOCUMENT
///
/// \relates region
bool cath::chop::has_chain_label(const region &prm_region ///< TODOCUMENT
                                  ) {
	return static_cast<bool>( prm_region.get_opt_chain_label() );
}

/// \brief TODOCUMENT
///
/// \relates region
const chain_label & cath::chop::get_chain_label(const region &prm_region ///< TODOCUMENT
                                                ) {
	return prm_region.get_opt_chain_label().value();
}

/// \brief TODOCUMENT
///
/// \relates region
const residue_name_opt & cath::chop::get_opt_start_name(const region &prm_region ///< TODOCUMENT
                                                        ) {
	return prm_region.get_start_residue().get_opt_residue_name();
}

/// \brief TODOCUMENT
///
/// \relates region
const residue_name_opt & cath::chop::get_opt_stop_name(const region &prm_region ///< TODOCUMENT
                                                       ) {
	return prm_region.get_stop_residue().get_opt_residue_name();
}

/// \brief TODOCUMENT
///
/// \relates region
const size_opt & cath::chop::get_opt_start_index(const region &prm_region ///< TODOCUMENT
                                                 ) {
	return prm_region.get_start_residue().get_opt_residue_index();
}

/// \brief TODOCUMENT
///
/// \relates region
const size_opt & cath::chop::get_opt_stop_index(const region &prm_region ///< TODOCUMENT
                                                ) {
	return prm_region.get_stop_residue().get_opt_residue_index();
}

/// \brief TODOCUMENT
///
/// \relates region
bool cath::chop::has_names(const region &prm_region ///< TODOCUMENT
                           ) {
	return prm_region.has_starts_stops() && static_cast<bool>( get_opt_start_name( prm_region ) );
}

/// \brief TODOCUMENT
///
/// \relates region
bool cath::chop::has_indices(const region &prm_region ///< TODOCUMENT
                             ) {
	return prm_region.has_starts_stops() && static_cast<bool>( get_opt_start_index( prm_region ) );
}

/// \brief TODOCUMENT
///
/// \relates region
void cath::chop::check_has_names(const region &prm_region ///< TODOCUMENT
                                 ) {
	if ( ! has_names ( prm_region ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to use an index-only region for calculations requiring residue names"));
	}
}

/// \brief TODOCUMENT
///
/// \relates region
void cath::chop::check_has_indices(const region &prm_region ///< TODOCUMENT
                                   ) {
	if ( ! has_indices( prm_region ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to use an name-only region for calculations requiring indices"));
	}
}

/// \brief TODOCUMENT
///
/// \relates region
const residue_name & cath::chop::get_start_name(const region &prm_region ///< TODOCUMENT
                                                ) {
	return *get_opt_start_name( prm_region );
}

/// \brief TODOCUMENT
///
/// \relates region
const residue_name & cath::chop::get_stop_name(const region &prm_region ///< TODOCUMENT
                                               ) {
	return *get_opt_stop_name( prm_region );
}

/// \brief TODOCUMENT
///
/// \relates region
const size_t & cath::chop::get_start_index(const region &prm_region ///< TODOCUMENT
                                           ) {
	return *get_opt_start_index( prm_region );
}

/// \brief TODOCUMENT
///
/// \relates region
const size_t & cath::chop::get_stop_index(const region &prm_region ///< TODOCUMENT
                                          ) {
	return *get_opt_stop_index( prm_region );
}

/// \brief Get the residue_id of the residue on which the specified region starts
///
/// \pre has_chain_label( prm_region ) && has_names( prm_region )
///
/// \relates region
///
/// \param prm_region The region to query
residue_id cath::chop::get_start_id( const region &prm_region ) {
	return { *prm_region.get_opt_chain_label(), get_start_name( prm_region ) };
}

/// \brief Get the residue_id of the residue on which the specified region stops
///
/// \pre has_chain_label( prm_region ) && has_names( prm_region )
///
/// \relates region
///
/// \param prm_region The region to query
residue_id cath::chop::get_stop_id( const region &prm_region ) {
	return { *prm_region.get_opt_chain_label(), get_stop_name( prm_region ) };
}

/// \brief TODOCUMENT
///
/// \relates region
residue_locating_opt cath::chop::get_residue_locating(const region &prm_region ///< TODOCUMENT
                                                      ) {
	return if_then_optional(
		prm_region.has_starts_stops(),
		get_residue_locating( prm_region.get_start_residue() )
	);
}

/// \brief TODOCUMENT
///
/// \relates region
size_t cath::chop::get_length(const region &prm_region ///< TODOCUMENT
                              ) {
	check_has_indices( prm_region );
	return 1 + get_stop_index ( prm_region )
	         - get_start_index( prm_region );
}

/// \brief TODOCUMENT
///
/// \relates region
region_comparison cath::chop::compare_locations(const region &prm_region_a, ///< TODOCUMENT
                                                const region &prm_region_b  ///< TODOCUMENT
                                                ) {
	// Check that both regions have indices
	check_has_indices( prm_region_a );
	check_has_indices( prm_region_b );

	// Grab the start and stop indices
	const size_t start_a = get_start_index( prm_region_a );
	const size_t start_b = get_start_index( prm_region_b );
	const size_t stop_a  = get_stop_index ( prm_region_a );
	const size_t stop_b  = get_stop_index ( prm_region_b );

	//      If the starts and stops are equal,                      then THE_SAME_AS
	if ( start_a == start_b && stop_a == stop_b ) { return region_comparison::THE_SAME_AS   ; }
	// Else if the first's start is no earlier and stop no later,   then A_SUBSET_OF
	if ( start_a >= start_b && stop_a <= stop_b ) { return region_comparison::A_SUBSET_OF   ; }
	// Else if the first's start is no later   and stop no earlier, then A_SUPERSET_OF
	if ( start_a <= start_b && stop_a >= stop_b ) { return region_comparison::A_SUPERSET_OF ; }
	// Otherwise, one of the regions is earlier than the other, even if overlappingly
	//
	// If   the first  is earlier then... If the first stops  before the second starts : STRICTLY_BEFORE else OVERLAPPINGLY_BEFORE
	if ( start_a < start_b ) { return ( stop_a  < start_b ) ? region_comparison::STRICTLY_BEFORE : region_comparison::OVERLAPPINGLY_BEFORE ; }
	// Else the second is earlier, so...  If the first starts after  the second stops  : STRICTLY_AFTER  else region_comparison::OVERLAPPINGLY_AFTER
	                           return ( start_a > stop_b  ) ? region_comparison::STRICTLY_AFTER  : region_comparison::OVERLAPPINGLY_AFTER  ;
}

/// \brief Expand the specified region to the whole chain of the region
///
/// `has_chain_label( prm_region )`
region cath::chop::expand_to_chain(const region &prm_region ///< The region to expand
                                   ) {
	if ( ! has_chain_label( prm_region ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot expand_to_chain() a region that doesn't have a chain label"));
	}
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return region{ get_chain_label( prm_region ) };
}

/// \brief Make a simple, whole-chain region from a chain label char
///
/// \relates region
region cath::chop::make_simple_region(const char &prm_chain_label_char ///< The chain label char for the new region
                                      ) {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return region{ chain_label { prm_chain_label_char } };
}

/// \brief Make a simple region from a start index and stop index
///
/// This is a direct pass-through to one of the region ctors and is just
/// present so that the make_simple_region() overload set covers all the
/// simple argument combinations.
///
/// \relates region
region cath::chop::make_simple_region(const size_t &prm_start_index, ///< The start index
                                      const size_t &prm_stop_index   ///< The stop index
                                      ) {
	return { prm_start_index, prm_stop_index };
}

/// \brief Make a simple region from a chain label char and a start & stop res name num
///
/// \relates region
region cath::chop::make_simple_region(const char &prm_chain_label_char,   ///< The chain label char for the new region
                                      const int  &prm_start_res_name_num, ///< The num for the start residue name
                                      const int  &prm_stop_res_name_num   ///< The num for the stop  residue name
                                      ) {
	return {
		chain_label { prm_chain_label_char   },
		residue_name{ prm_start_res_name_num },
		residue_name{ prm_stop_res_name_num  }
	};
}

/// \brief Make a simple region from a chain label char and a start & stop res name num with insert code chars
///
/// \relates region
region cath::chop::make_simple_region(const char &prm_chain_label_char,   ///< The chain label char for the new region
                                      const int  &prm_start_res_name_num, ///< The num for the start residue name
                                      const char &prm_start_res_name_ins, ///< The insert code char for the start  residue name
                                      const int  &prm_stop_res_name_num,  ///< The num for the stop  residue name
                                      const char &prm_stop_res_name_ins   ///< The insert code char for the stop  residue name
                                      ) {
	return {
		chain_label { prm_chain_label_char   },
		residue_name{ prm_start_res_name_num, prm_start_res_name_ins },
		residue_name{ prm_stop_res_name_num,  prm_stop_res_name_ins  }
	};
}

/// \brief Generate a string describing the specified region
///
/// \relates region
string cath::chop::to_string(const region &prm_region ///< The region to describe
                             ) {
	str_vec parts;
	if ( has_chain_label( prm_region ) ) {
		parts.push_back( "chain:"      +      get_chain_label( prm_region ).to_string()  );
	}
	if ( has_names( prm_region ) ) {
		parts.push_back( "start_name:" +      to_string( get_start_name ( prm_region ) ) );
		parts.push_back( "stop_name:"  +      to_string( get_stop_name  ( prm_region ) ) );
	}
	if ( has_indices( prm_region ) ) {
		parts.push_back( "start_idx:"  + std::to_string( get_start_index ( prm_region ) ) );
		parts.push_back( "stop_idx:"   + std::to_string( get_stop_index  ( prm_region ) ) );
	}
	return "region[ " + join( parts, ", " ) + " ]";
}

/// \brief Insert a description of the specified region into the specified ostream
///
/// \relates region
ostream & cath::chop::operator<<(ostream      &prm_os,    ///< The ostream into which the description should be inserted
                                 const region &prm_region ///< The region to describe
                                 ) {
	prm_os << to_string( prm_region );
	return prm_os;
}

/// \brief Generate a string describing the specified region_vec
///
/// \relates region
string cath::chop::to_string(const region_vec &prm_regions ///< The region_vec to describe
                             ) {
	return join(
		prm_regions
			| transformed( [] (const region &x) { return to_string( x ); } ),
		","
	);
}

/// \brief Insert a description of the specified region_vec into the specified ostream
///
/// \relates region_vec
ostream & cath::chop::operator<<(ostream      &prm_os,    ///< The ostream into which the description should be inserted
                                 const region_vec &prm_regions ///< The region_vec to describe
                                 ) {
	prm_os << to_string( prm_regions );
	return prm_os;
}

/// \brief TODOCUMENT
void cath::chop::validate(any           &prm_prev_value,    ///< TODOCUMENT
                          const str_vec &prm_value_strings, ///< TODOCUMENT
                          domain *,
                          int
                          ) {
	prm_prev_value = [&] {
		// Standard validate boilerplate:
		//  * Make sure no previous assignment to 'a' was made.
		//  * Extract the first string from 'prm_value_strings'.
		//    (If there is more than one string, it's an error, and exception will be thrown.)
		boost::program_options::validators::check_first_occurrence( prm_prev_value );
		const std::string &value_string = boost::program_options::validators::get_single_string( prm_value_strings );

		// Attempt to lexical_cast value_string and if it fails, throw an invalid_option_value exception
		try {
			return sillitoe_chopping_format{}.parse_domain( value_string );
		}
		catch (...) {
			BOOST_THROW_EXCEPTION( boost::program_options::invalid_option_value( value_string ) );
		}
	} ();
}
