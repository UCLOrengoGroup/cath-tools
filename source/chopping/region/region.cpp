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

#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::chop;
using namespace cath::common;

using boost::algorithm::join;
using std::ostream;
using std::string;
using std::make_pair;

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
region::region(const chain_label &arg_chain_label ///< The chain_label for the whole-chain region to construct
               ) : the_chain_label( arg_chain_label ) {
}

/// \brief Ctor for region
region::region(const chain_label  &arg_chain_label,        ///<TODOCUMENT
               const residue_name &arg_start_residue_name, ///<TODOCUMENT
               const residue_name &arg_stop_residue_name   ///<TODOCUMENT
               ) : the_chain_label( arg_chain_label        ),
                   residues       { make_pair(
                   	residue_location{ arg_start_residue_name },
                   	residue_location{ arg_stop_residue_name  }
                   ) } {
	sanity_check();
}

/// \brief Ctor for region
region::region(const chain_label  &arg_chain_label,        ///<TODOCUMENT
               const residue_name &arg_start_residue_name, ///<TODOCUMENT
               const size_t       &arg_start_index,        ///<TODOCUMENT
               const residue_name &arg_stop_residue_name,  ///<TODOCUMENT
               const size_t       &arg_stop_index          ///<TODOCUMENT
               ) : the_chain_label( arg_chain_label                         ),
                   residues       { make_pair(
                   	residue_location{ arg_start_residue_name, arg_start_index },
                   	residue_location{ arg_stop_residue_name,  arg_stop_index  }
                   ) } {
	sanity_check();
}

/// \brief Ctor for region
region::region(const size_t &arg_start_index, ///<TODOCUMENT
               const size_t &arg_stop_index   ///<TODOCUMENT
               ) : residues{ make_pair(
                   	residue_location{ arg_start_index },
                   	residue_location{ arg_stop_index  }
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
bool cath::chop::operator==(const region &arg_lhs, ///< The first  region to compare
                            const region &arg_rhs  ///< The second region to compare
                            ) {
	return (
		( arg_lhs.get_opt_chain_label() == arg_rhs.get_opt_chain_label() )
		&&
		( arg_lhs.has_starts_stops()    == arg_rhs.has_starts_stops() )
		&&
		(
			arg_lhs.has_starts_stops()
				? (
					arg_lhs.get_start_residue()   == arg_rhs.get_start_residue()
					&&
					arg_lhs.get_stop_residue()    == arg_rhs.get_stop_residue()
				)
				: true
		)
	);
}

/// \brief TODOCUMENT
///
/// \relates region
bool cath::chop::has_chain_label(const region &arg_region ///< TODOCUMENT
                                  ) {
	return static_cast<bool>( arg_region.get_opt_chain_label() );
}

/// \brief TODOCUMENT
///
/// \relates region
const chain_label & cath::chop::get_chain_label(const region &arg_region ///< TODOCUMENT
                                                ) {
	return arg_region.get_opt_chain_label().get();
}

/// \brief TODOCUMENT
///
/// \relates region
const residue_name_opt & cath::chop::get_opt_start_name(const region &arg_region ///< TODOCUMENT
                                                        ) {
	return arg_region.get_start_residue().get_opt_residue_name();
}

/// \brief TODOCUMENT
///
/// \relates region
const residue_name_opt & cath::chop::get_opt_stop_name(const region &arg_region ///< TODOCUMENT
                                                       ) {
	return arg_region.get_stop_residue().get_opt_residue_name();
}

/// \brief TODOCUMENT
///
/// \relates region
const size_opt & cath::chop::get_opt_start_index(const region &arg_region ///< TODOCUMENT
                                                 ) {
	return arg_region.get_start_residue().get_opt_residue_index();
}

/// \brief TODOCUMENT
///
/// \relates region
const size_opt & cath::chop::get_opt_stop_index(const region &arg_region ///< TODOCUMENT
                                                ) {
	return arg_region.get_stop_residue().get_opt_residue_index();
}

/// \brief TODOCUMENT
///
/// \relates region
bool cath::chop::has_names(const region &arg_region ///< TODOCUMENT
                           ) {
	return arg_region.has_starts_stops() && static_cast<bool>( get_opt_start_name( arg_region ) );
}

/// \brief TODOCUMENT
///
/// \relates region
bool cath::chop::has_indices(const region &arg_region ///< TODOCUMENT
                             ) {
	return arg_region.has_starts_stops() && static_cast<bool>( get_opt_start_index( arg_region ) );
}

/// \brief TODOCUMENT
///
/// \relates region
void cath::chop::check_has_names(const region &arg_region ///< TODOCUMENT
                                 ) {
	if ( ! has_names ( arg_region ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to use an index-only region for calculations requiring residue names"));
	}
}

/// \brief TODOCUMENT
///
/// \relates region
void cath::chop::check_has_indices(const region &arg_region ///< TODOCUMENT
                                   ) {
	if ( ! has_indices( arg_region ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to use an name-only region for calculations requiring indices"));
	}
}

/// \brief TODOCUMENT
///
/// \relates region
const residue_name & cath::chop::get_start_name(const region &arg_region ///< TODOCUMENT
                                                ) {
	return *get_opt_start_name( arg_region );
}

/// \brief TODOCUMENT
///
/// \relates region
const residue_name & cath::chop::get_stop_name(const region &arg_region ///< TODOCUMENT
                                               ) {
	return *get_opt_stop_name( arg_region );
}

/// \brief TODOCUMENT
///
/// \relates region
const size_t & cath::chop::get_start_index(const region &arg_region ///< TODOCUMENT
                                           ) {
	return *get_opt_start_index( arg_region );
}

/// \brief TODOCUMENT
///
/// \relates region
const size_t & cath::chop::get_stop_index(const region &arg_region ///< TODOCUMENT
                                          ) {
	return *get_opt_stop_index( arg_region );
}

/// \brief TODOCUMENT
///
/// \relates region
residue_locating cath::chop::get_residue_locating(const region &arg_region ///< TODOCUMENT
                                      ) {
	return get_residue_locating( arg_region.get_start_residue() );
}

/// \brief TODOCUMENT
///
/// \relates region
size_t cath::chop::get_length(const region &arg_region ///< TODOCUMENT
                              ) {
	check_has_indices( arg_region );
	return 1 + get_stop_index ( arg_region )
	         - get_start_index( arg_region );
}

/// \brief TODOCUMENT
///
/// \relates region
region_comparison cath::chop::compare_locations(const region &arg_region_a, ///< TODOCUMENT
                                                const region &arg_region_b  ///< TODOCUMENT
                                                ) {
	// Check that both regions have indices
	check_has_indices( arg_region_a );
	check_has_indices( arg_region_b );

	// Grab the start and stop indices
	const size_t start_a = get_start_index( arg_region_a );
	const size_t start_b = get_start_index( arg_region_b );
	const size_t stop_a  = get_stop_index ( arg_region_a );
	const size_t stop_b  = get_stop_index ( arg_region_b );

	//      If the starts and stops are equal,                      then THE_SAME_AS
	if      ( start_a == start_b && stop_a == stop_b ) { return region_comparison::THE_SAME_AS   ; }
	// Else if the first's start is no earlier and stop no later,   then A_SUBSET_OF
	else if ( start_a >= start_b && stop_a <= stop_b ) { return region_comparison::A_SUBSET_OF   ; }
	// Else if the first's start is no later   and stop no earlier, then A_SUPERSET_OF
	else if ( start_a <= start_b && stop_a >= stop_b ) { return region_comparison::A_SUPERSET_OF ; }
	// Otherwise, one of the regions is earlier than the other, even if overlappingly
	//
	// If   the first  is earlier then... If the first stops  before the second starts : STRICTLY_BEFORE else OVERLAPPINGLY_BEFORE
	else if ( start_a < start_b ) { return ( stop_a  < start_b ) ? region_comparison::STRICTLY_BEFORE : region_comparison::OVERLAPPINGLY_BEFORE ; }
	// Else the second is earlier, so...  If the first starts after  the second stops  : STRICTLY_AFTER  else region_comparison::OVERLAPPINGLY_AFTER
	else                          { return ( start_a > stop_b  ) ? region_comparison::STRICTLY_AFTER  : region_comparison::OVERLAPPINGLY_AFTER  ; }
}

/// \brief Make a simple, whole-chain region from a chain label char
///
/// \relates region
region cath::chop::make_simple_region(const char &arg_chain_label_char ///< The chain label char for the new region
                                      ) {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return region{ chain_label { arg_chain_label_char } };
}

/// \brief Make a simple region from a start index and stop index
///
/// This is a direct pass-through to one of the region ctors and is just
/// present so that the make_simple_region() overload set covers all the
/// simple argument combinations.
///
/// \relates region
region cath::chop::make_simple_region(const size_t &arg_start_index, ///< The start index
                                      const size_t &arg_stop_index   ///< The stop index
                                      ) {
	return { arg_start_index, arg_stop_index };
}

/// \brief Make a simple region from a chain label char and a start & stop res name num
///
/// \relates region
region cath::chop::make_simple_region(const char &arg_chain_label_char,   ///< The chain label char for the new region
                                      const int  &arg_start_res_name_num, ///< The num for the start residue name
                                      const int  &arg_stop_res_name_num   ///< The num for the stop  residue name
                                      ) {
	return {
		chain_label { arg_chain_label_char   },
		residue_name{ arg_start_res_name_num },
		residue_name{ arg_stop_res_name_num  }
	};
}

/// \brief Make a simple region from a chain label char and a start & stop res name num with insert code chars
///
/// \relates region
region cath::chop::make_simple_region(const char &arg_chain_label_char,   ///< The chain label char for the new region
                                      const int  &arg_start_res_name_num, ///< The num for the start residue name
                                      const char &arg_start_res_name_ins, ///< The insert code char for the start  residue name
                                      const int  &arg_stop_res_name_num,  ///< The num for the stop  residue name
                                      const char &arg_stop_res_name_ins   ///< The insert code char for the stop  residue name
                                      ) {
	return {
		chain_label { arg_chain_label_char   },
		residue_name{ arg_start_res_name_num, arg_start_res_name_ins },
		residue_name{ arg_stop_res_name_num,  arg_stop_res_name_ins  }
	};
}

/// \brief Generate a string describing the specified region
///
/// \relates region
string cath::chop::to_string(const region &arg_region ///< The region to describe
                             ) {
	str_vec parts;
	if ( has_chain_label( arg_region ) ) {
		parts.push_back( "chain:"      +      to_string( get_chain_label( arg_region ) ) );
	}
	if ( has_names( arg_region ) ) {
		parts.push_back( "start_name:" +      to_string( get_start_name ( arg_region ) ) );
		parts.push_back( "stop_name:"  +      to_string( get_stop_name  ( arg_region ) ) );
	}
	if ( has_indices( arg_region ) ) {
		parts.push_back( "start_idx:"  + std::to_string( get_start_index ( arg_region ) ) );
		parts.push_back( "stop_idx:"   + std::to_string( get_stop_index  ( arg_region ) ) );
	}
	return "region{ " + join( parts, ", " ) + " }";
}

/// \brief Insert a description of the specified region into the specified ostream
///
/// \relates region
ostream & cath::chop::operator<<(ostream      &arg_os,    ///< The ostream into which the description should be inserted
                                 const region &arg_region ///< The region to describe
                                 ) {
	arg_os << to_string( arg_region );
	return arg_os;
}
