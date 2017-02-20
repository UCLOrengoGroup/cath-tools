/// \file
/// \brief The dssp_ss_calc class definitions

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

#include "dssp_ss_calc.hpp"

#include <boost/algorithm/clamp.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/append.hpp"
#include "common/algorithm/transform_build.hpp"
// #include "common/boost_addenda/range/adaptor/lexical_casted.hpp" // ***** TEMPORARY *****
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc_type.hpp"
#include "structure/sec_struc_calc/dssp/bifur_hbond_list.hpp"
#include "structure/sec_struc_calc/dssp/dssp_hbond_calc.hpp"

#include <algorithm>

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::sec;
using namespace cath::sec::detail;

using boost::adaptors::reversed;
using boost::algorithm::any_of;
using boost::algorithm::clamp;
using boost::algorithm::join;
using boost::format;
using boost::irange;
using boost::make_optional;
using boost::none;
using boost::optional;
using std::max;
using std::min;
using std::ostream;
using std::string;
using std::vector;

constexpr size_t sec_struc_consts::MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE;
constexpr size_t sec_struc_consts::BETA_BULGE_MAX_DIFF_SOURCE;
constexpr size_t sec_struc_consts::BETA_BULGE_MAX_DIFF_DEST;
constexpr beta_bridge_context beta_bridge::DEFAULT_CONTEXT;

/// \brief Generate a string describing the specified beta_bridge_type
///
/// \relates beta_bridge_type
string cath::sec::detail::to_string(const beta_bridge_type &arg_beta_bridge_type ///< The beta_bridge_type to describe
                                    ) {
	switch ( arg_beta_bridge_type ) {
		case ( beta_bridge_type::PARALLEL      ) : { return "PARALLEL"; };
		case ( beta_bridge_type::ANTI_PARALLEL ) : { return "ANTIPARA"; };
		default : {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of beta_bridge_type not recognised whilst converting to_string()"));
			return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}
}

/// \brief Insert a description of the specified beta_bridge_type into the specified ostream
///
/// \relates beta_bridge_type
ostream & cath::sec::detail::operator<<(ostream                &arg_os,              ///< The ostream into which the description should be inserted
                                        const beta_bridge_type &arg_beta_bridge_type ///< The beta_bridge_type to describe
                                        ) {
	arg_os << to_string( arg_beta_bridge_type );
	return arg_os;
}

/// \brief Generate a string describing the specified beta_bridge_context
///
/// \relates beta_bridge_context
string cath::sec::detail::to_string(const beta_bridge_context &arg_beta_bridge_context ///< The beta_bridge_context to describe
                                    ) {
	switch ( arg_beta_bridge_context ) {
		case ( beta_bridge_context::LONE_BRIDGE ) : { return "BRDGE"; };
		case ( beta_bridge_context::IN_SHEET    ) : { return "SHEET"; };
		default : {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of beta_bridge_context not recognised whilst converting to_string()"));
			return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}
}

/// \brief Insert a description of the specified beta_bridge_context into the specified ostream
///
/// \relates beta_bridge_context
ostream & cath::sec::detail::operator<<(ostream                   &arg_os,                 ///< The ostream into which the description should be inserted
                                        const beta_bridge_context &arg_beta_bridge_context ///< The beta_bridge_context to describe
                                        ) {
	arg_os << to_string( arg_beta_bridge_context );
	return arg_os;
}

/// \brief Return whether the two specified beta_bridges are identical
///
/// \relates beta_bridge
inline bool cath::sec::detail::operator==(const beta_bridge &arg_lhs, ///< The first  beta_bridge to compare
                                          const beta_bridge &arg_rhs  ///< The second beta_bridge to compare
                                          ) {
	return (
		arg_lhs.partner_idx == arg_rhs.partner_idx
		&&
		arg_lhs.type        == arg_rhs.type
	);
}

/// \brief Generate a string describing the specified beta_bridge
///
/// \relates beta_bridge
string cath::sec::detail::to_string(const beta_bridge &arg_beta_bridge ///< The beta_bridge to describe
                                    ) {
	return "beta_bridge[partner:"
		+ ( format( "%5g") % arg_beta_bridge.partner_idx ).str()
		+ ", type:"
		+ to_string( arg_beta_bridge.type )
		+ ", context:"
		+ to_string( arg_beta_bridge.context )
		+ "]";
}

/// \brief Insert a description of the specified beta_bridge into the specified ostream
///
/// \relates beta_bridge
ostream & cath::sec::detail::operator<<(ostream           &arg_os,         ///< The ostream into which the description should be inserted
                                        const beta_bridge &arg_beta_bridge ///< The beta_bridge to describe
                                        ) {
	arg_os << to_string( arg_beta_bridge );
	return arg_os;
}

/// \brief Return whether the specified helix_category value means that the residue in question is helix-bonded to an earlier residue
///
/// \relates helix_category
bool cath::sec::detail::is_bonded_to_earlier(const helix_category &arg_helix_cat ///<  The helix_category to query
                                             ) {
	return ( arg_helix_cat == helix_category::BONDED_TO_BOTH || arg_helix_cat == helix_category::BONDED_TO_EARLIER_ONLY );
}

/// \brief Return whether the specified helix_category value means that the residue in question is helix-bonded to a later residue
///
/// \relates helix_category
bool cath::sec::detail::is_bonded_to_later(const helix_category &arg_helix_cat ///< The helix_category to query
                                           ) {
	return ( arg_helix_cat == helix_category::BONDED_TO_BOTH || arg_helix_cat == helix_category::BONDED_TO_LATER_ONLY );
}


/// \brief Return whether the specified hbond_half_opt_pair is bonded to the residue at the specified index
bool cath::sec::detail::is_bonded_to(const hbond_half_opt_pair &arg_bound_pair, ///< The hbond_half_opt_pair to examine
                                     const size_t              &arg_dest_index  ///< The index of the residue to which we're interested in potential bonds
                                     ) {
	return (
		is_bondy_enough( arg_bound_pair.first )
		&&
		(
			( arg_bound_pair.first->index == arg_dest_index )
			||
			(
				is_bondy_enough( arg_bound_pair.second )
				&&
				arg_bound_pair.second->index == arg_dest_index
			)
		)
	);
}

/// \brief Get the n-helix category of the specified bifur_hbond at the specified index
optional<helix_category> cath::sec::detail::n_helix_cat(const bifur_hbond &arg_bifur_hbond, ///< The bifur_hbond to query
                                                        const size_t      &arg_index,       ///< The index of the bifur_hbond
                                                        const size_t      &arg_helix_num    ///< "n", the difference between bonded residues
                                                        ) {
	const bool bonded_to_later   = is_bonded_to( arg_bifur_hbond.get_bound_pair_for_this_co(), arg_index + arg_helix_num );
	const bool bonded_to_earlier = ( arg_index >= arg_helix_num
	                               &&
	                               is_bonded_to( arg_bifur_hbond.get_bound_pair_for_this_nh(), arg_index - arg_helix_num ) );

	if ( bonded_to_later ) {
		return { bonded_to_earlier ? helix_category::BONDED_TO_BOTH : helix_category::BONDED_TO_LATER_ONLY };
	}
	else {
		return make_optional( bonded_to_earlier, helix_category::BONDED_TO_EARLIER_ONLY );
	}
}

/// \brief Get the four-helix category of the bifur_hbond at the specified index of the specified bifur_hbond_list
optional<helix_category> cath::sec::detail::n_helix_cat(const bifur_hbond_list &arg_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                        const size_t           &arg_index,            ///< The index of the bifur_hbond in question
                                                        const size_t           &arg_helix_num         ///< "n", the difference between bonded residues
                                                        ) {
	return n_helix_cat( arg_bifur_hbond_list[ arg_index ], arg_index, arg_helix_num );
}

/// \brief Return whether the residue at the specified index is within a four-helix according to the specified bifur_hbond_list
bool cath::sec::detail::is_n_helix(const bifur_hbond_list &arg_bifur_hbond_list, ///< The bifur_hbond_list to query
                                   const size_t           &arg_index,            ///< The index of the residue in question
                                   const size_t           &arg_helix_num         ///< "n", the difference between bonded residues
                                   ) {
	for (const size_t &second_index : irange( max( arg_helix_num, arg_index ) - ( arg_helix_num - 1_z ), arg_index + 1_z ) | reversed ) {
		const auto first_cat  = n_helix_cat( arg_bifur_hbond_list, second_index - 1, arg_helix_num );
		const auto second_cat = n_helix_cat( arg_bifur_hbond_list, second_index,     arg_helix_num );
		if ( first_cat && is_bonded_to_later( *first_cat ) && second_cat && is_bonded_to_later( *second_cat ) ) {
			return true;
		}
	}
	return false;
	// // Grab the helix_category of the specified index in the specified bifur_hbond_list
	// const auto fhc = n_helix_cat( arg_bifur_hbond_list, arg_index );
	// if ( fhc ) {
	// 	// If bonded to both, then true
	// 	if ( *fhc == helix_category::BONDED_TO_BOTH ) {
	// 		return true;
	// 	}
	// 	// If bonded to later then true if followed by another that's bonded to later
	// 	else if ( *fhc == helix_category::BONDED_TO_LATER_ONLY && arg_index > 0 ) {
	// 		const auto prev_cat = n_helix_cat( arg_bifur_hbond_list, arg_index - 1 );
	// 		if ( prev_cat && ( *prev_cat == *fhc || *prev_cat == helix_category::BONDED_TO_BOTH ) ) {
	// 			return true;
	// 		}
	// 	}
	// 	// If bonded to earlier then true if preceded by another that's bonded to earlier
	// 	else if ( *fhc == helix_category::BONDED_TO_EARLIER_ONLY && arg_index + 1 < arg_bifur_hbond_list.size() ) {
	// 		const auto next_cat = n_helix_cat( arg_bifur_hbond_list, arg_index + 1 );
	// 		if ( next_cat && ( *next_cat == *fhc || *next_cat == helix_category::BONDED_TO_BOTH ) ) {
	// 			return true;
	// 		}
	// 	}
	// }
	// return false;
}

/// \brief Return whether the specified index is within range to be the middle of a beta triplet
bool cath::sec::detail::beta_index_in_range(const bifur_hbond_list &arg_bifur_hbond_list, ///< The bifur_hbond_list providing the context
                                            const size_t           &arg_index             ///< The index to check
                                            ) {
	return (
		( arg_index > 0 )
		&&
		( arg_index + 1 < arg_bifur_hbond_list.size() )
	);
}

/// \brief Return whether there are parallel-beta-bridge bonds between the specified indices
///        that bond to the source residue itself
beta_bridge_opt cath::sec::detail::has_parallel_beta_bridge_bonds_to_src(const bifur_hbond_list &arg_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                         const size_t           &arg_src_index,        ///< The index of the source residue to query
                                                                         const size_t           &arg_dest_index        ///< The index of the destination residue to query
                                                                         ) {
	return make_optional(
		(
			difference( arg_src_index, arg_dest_index ) >= sec_struc_consts::MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE
			&&
			beta_index_in_range( arg_bifur_hbond_list, arg_src_index  )
			&&
			beta_index_in_range( arg_bifur_hbond_list, arg_dest_index )
			&&
			is_bonded_to( arg_bifur_hbond_list[ arg_src_index     ].get_bound_pair_for_this_nh(), arg_dest_index - 1 )
			&&
			is_bonded_to( arg_bifur_hbond_list[ arg_src_index     ].get_bound_pair_for_this_co(), arg_dest_index + 1 )
		),
		beta_bridge{ arg_dest_index, beta_bridge_type::PARALLEL }
	);
}

/// \brief Return whether there are parallel-beta-bridge bonds between the specified indices
///        that bond to the residues straddling the source residue
beta_bridge_opt cath::sec::detail::has_parallel_beta_bridge_bonds_straddling_src(const bifur_hbond_list &arg_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                                 const size_t           &arg_src_index,        ///< The index of the source residue to query
                                                                                 const size_t           &arg_dest_index        ///< The index of the destination residue to query
                                                                                 ) {
	return make_optional(
		(
			difference( arg_src_index, arg_dest_index ) >= sec_struc_consts::MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE
			&&
			beta_index_in_range( arg_bifur_hbond_list, arg_src_index  )
			&&
			beta_index_in_range( arg_bifur_hbond_list, arg_dest_index )
			&&
			is_bonded_to( arg_bifur_hbond_list[ arg_src_index - 1 ].get_bound_pair_for_this_co(), arg_dest_index    )
			&&
			is_bonded_to( arg_bifur_hbond_list[ arg_src_index + 1 ].get_bound_pair_for_this_nh(), arg_dest_index    )
		),
		beta_bridge{ arg_dest_index, beta_bridge_type::PARALLEL }
	);
}

/// \brief Return whether there are antiparallel-beta-bridge bonds between the specified indices
///        that bond to the source residue itself
beta_bridge_opt cath::sec::detail::has_antiparallel_beta_bridge_bonds_to_src(const bifur_hbond_list &arg_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                             const size_t           &arg_src_index,        ///< The index of the source residue to query
                                                                             const size_t           &arg_dest_index        ///< The index of the destination residue to query
                                                                             ) {
	return make_optional(
		(
			difference( arg_src_index, arg_dest_index ) >= sec_struc_consts::MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE
			&&
			arg_src_index  < arg_bifur_hbond_list.size()
			&&
			arg_dest_index < arg_bifur_hbond_list.size()
			&&
			is_bonded_to( arg_bifur_hbond_list[ arg_src_index     ].get_bound_pair_for_this_nh(), arg_dest_index     )
			&&
			is_bonded_to( arg_bifur_hbond_list[ arg_src_index     ].get_bound_pair_for_this_co(), arg_dest_index     )
		),
		beta_bridge{ arg_dest_index, beta_bridge_type::ANTI_PARALLEL }
	);
}

/// \brief Return whether there are antiparallel-beta-bridge bonds between the specified indices
///        that bond to the residues straddling the source residue
beta_bridge_opt cath::sec::detail::has_antiparallel_beta_bridge_bonds_straddling_src(const bifur_hbond_list &arg_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                                     const size_t           &arg_src_index,        ///< The index of the source residue to query
                                                                                     const size_t           &arg_dest_index        ///< The index of the destination residue to query
                                                                                     ) {
	return make_optional(
		(
			difference( arg_src_index, arg_dest_index ) >= sec_struc_consts::MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE
			&&
			beta_index_in_range( arg_bifur_hbond_list, arg_src_index  )
			&&
			beta_index_in_range( arg_bifur_hbond_list, arg_dest_index )
			&&
			is_bonded_to( arg_bifur_hbond_list[ arg_src_index - 1 ].get_bound_pair_for_this_co(), arg_dest_index + 1 )
			&&
			is_bonded_to( arg_bifur_hbond_list[ arg_src_index + 1 ].get_bound_pair_for_this_nh(), arg_dest_index - 1 )
		),
		beta_bridge{ arg_dest_index, beta_bridge_type::ANTI_PARALLEL }
	);
}

/// \brief Return whether the residue at the specified index is part of a parallel beta bridge according to the specified bifur_hbond_list
beta_bridge_vec cath::sec::detail::has_parallel_beta_bridge(const bifur_hbond_list &arg_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                            const size_t           &arg_index             ///< The index to query
                                                            ) {
	if ( ! beta_index_in_range( arg_bifur_hbond_list, arg_index ) ) {
		return {};
	}

	beta_bridge_vec results;
	const auto cnsdr_res_fn = [&] (const beta_bridge_opt &x) {
		if ( x ) {
			if ( ! contains( results, *x ) ) {
				results.push_back( *x );
			}
		}
	};
	const auto to_src_fn    = [&] (const size_t &x) { cnsdr_res_fn( has_parallel_beta_bridge_bonds_to_src        ( arg_bifur_hbond_list, arg_index, x ) ); };
	const auto strdl_src_fn = [&] (const size_t &x) { cnsdr_res_fn( has_parallel_beta_bridge_bonds_straddling_src( arg_bifur_hbond_list, arg_index, x ) ); };

	const auto prev_co_fn   = [&] (const auto   &x) { if ( is_bondy_enough( x )                     ) { strdl_src_fn( x->index     ); } };
	const auto this_co_fn   = [&] (const auto   &x) { if ( is_bondy_enough( x ) && ( x->index > 1 ) ) { to_src_fn   ( x->index - 1 ); } };
	const auto next_nh_fn   = [&] (const auto   &x) { if ( is_bondy_enough( x )                     ) { strdl_src_fn( x->index     ); } };

	const auto &prev_co = arg_bifur_hbond_list[ arg_index - 1 ].get_bound_pair_for_this_co();
	const auto &this_co = arg_bifur_hbond_list[ arg_index     ].get_bound_pair_for_this_co();
	const auto &next_nh = arg_bifur_hbond_list[ arg_index + 1 ].get_bound_pair_for_this_nh();

	prev_co_fn( prev_co.first  );
	prev_co_fn( prev_co.second );
	this_co_fn( this_co.first  );
	this_co_fn( this_co.second );
	next_nh_fn( next_nh.first  );
	next_nh_fn( next_nh.second );
	return results;
}

/// \brief Return whether the residue at the specified index is part of an antiparallel beta bridge according to the specified bifur_hbond_list
beta_bridge_vec cath::sec::detail::has_antiparallel_beta_bridge(const bifur_hbond_list &arg_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                const size_t           &arg_index             ///< The index to query
                                                                ) {
	if ( ! beta_index_in_range( arg_bifur_hbond_list, arg_index ) ) {
		return {};
	}

	beta_bridge_vec results;
	const auto cnsdr_res_fn = [&] (const beta_bridge_opt &x) {
		if ( x ) {
			if ( ! contains( results, *x ) ) {
				results.push_back( *x );
			}
		}
	};

	const auto sz = arg_bifur_hbond_list.size();

	const auto to_src_fn    = [&] (const size_t &x) { cnsdr_res_fn( has_antiparallel_beta_bridge_bonds_to_src        ( arg_bifur_hbond_list, arg_index, x ) ); };
	const auto strdl_src_fn = [&] (const size_t &x) { cnsdr_res_fn( has_antiparallel_beta_bridge_bonds_straddling_src( arg_bifur_hbond_list, arg_index, x ) ); };

	/// \todo Come GCCs all > 5.2.1, use a generic lambda here (causes ICE in GCC 5.2.1)
	const auto prev_co_fn   = [&] (const hbond_half_opt &x) { if ( is_bondy_enough( x ) && (   x->index       > 1  ) ) { strdl_src_fn( x->index - 1 ); } };
	const auto this_co_fn   = [&] (const hbond_half_opt &x) { if ( is_bondy_enough( x )                              ) { to_src_fn   ( x->index     ); } };
	const auto next_nh_fn   = [&] (const hbond_half_opt &x) { if ( is_bondy_enough( x ) && ( ( x->index + 1 ) < sz ) ) { strdl_src_fn( x->index + 1 ); } };

	const auto &prev_co = arg_bifur_hbond_list[ arg_index - 1 ].get_bound_pair_for_this_co();
	const auto &this_co = arg_bifur_hbond_list[ arg_index     ].get_bound_pair_for_this_co();
	const auto &next_nh = arg_bifur_hbond_list[ arg_index + 1 ].get_bound_pair_for_this_nh();

	prev_co_fn( prev_co.first  );
	prev_co_fn( prev_co.second );
	this_co_fn( this_co.first  );
	this_co_fn( this_co.second );
	next_nh_fn( next_nh.first  );
	next_nh_fn( next_nh.second );
	return results;
}

/// \brief Return whether the specified index is part of a beta bridge according to the specified bifur_hbond_list
beta_bridge_vec cath::sec::detail::has_beta_bridge(const bifur_hbond_list &arg_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                   const size_t           &arg_index             ///< The index to query
                                                   ) {
	return append_copy(
		has_parallel_beta_bridge    ( arg_bifur_hbond_list, arg_index ),
		has_antiparallel_beta_bridge( arg_bifur_hbond_list, arg_index )
	);
}

/// \brief Return whether there is a bulge between the residues at the other side of the specified beta bridges at the specified indices
///
/// A beta-bulge is where there are a few residues inserted in one strand of a beta-sheet, making it bulge out.
/// The residues in the bulge are still labelled as part of the strand.
///
/// The beta bridges must both be in the same direction for this to detect the bulge
///
/// \relates beta_bridge
bool cath::sec::detail::is_beta_bulge(const beta_bridge &arg_bridge_1, ///< The first  beta-bridge
                                      const size_t      &arg_index_1,  ///< The source index of the first  beta-bridge
                                      const beta_bridge &arg_bridge_2, ///< The second beta-bridge
                                      const size_t      &arg_index_2   ///< The source index of the second beta-bridge
                                      ) {
	const size_t src_index_diff  = difference( arg_index_1,              arg_index_2              );
	const size_t dest_index_diff = difference( arg_bridge_1.partner_idx, arg_bridge_2.partner_idx );
	const bool   src_ascending   = ( arg_index_2              > arg_index_1              );
	const bool   dest_ascending  = ( arg_bridge_2.partner_idx > arg_bridge_1.partner_idx );
	const bool   dirns_match     = ( src_ascending == dest_ascending );
	return (
		arg_bridge_1.type == arg_bridge_2.type
		&&
		clamp( src_index_diff,  1, sec_struc_consts::BETA_BULGE_MAX_DIFF_SOURCE ) == src_index_diff
		&&
		clamp( dest_index_diff, 1, sec_struc_consts::BETA_BULGE_MAX_DIFF_DEST   ) == dest_index_diff
		&&
		src_index_diff != dest_index_diff
		&&
		(
			( arg_bridge_1.type == beta_bridge_type::PARALLEL      &&   dirns_match )
			||
			( arg_bridge_2.type == beta_bridge_type::ANTI_PARALLEL && ! dirns_match )
		)
	);
}

/// \relates Add a beta-bulge between the specified residues in the specified sec_struc_types
///
/// Note: the from/to residues needn't necessarily be in ascending order
void cath::sec::detail::add_beta_bulge(sec_struc_type_vec &arg_secs, ///< The vector of sec_struc_type values to update
                                       const size_t       &arg_from, ///< The index of the residue at one end of the bulge
                                       const size_t       &arg_to    ///< The index of the residue at the other end of the bulge
                                       ) {
	for (const size_t &x : irange( min( arg_from, arg_to ), max( arg_from, arg_to ) + 1_z ) ) {
		arg_secs[ x ] = sec_struc_type::BETA_STRAND;
	}
}

/// \brief Return whether the two specified beta bridges are consecutive bridges in a beta sheet
///        (assuming the first's residue immediately follows the second's)
bool cath::sec::detail::are_consecutive_bridges_in_sheet(const beta_bridge &arg_beta_bridge_1, ///< The first  bridge to query
                                                         const beta_bridge &arg_beta_bridge_2  ///< The second bridge to query
                                                         ) {
	return (
		( arg_beta_bridge_1.type == arg_beta_bridge_2.type )
		&&
		( arg_beta_bridge_1.type != beta_bridge_type::PARALLEL      || arg_beta_bridge_1.partner_idx < arg_beta_bridge_2.partner_idx )
		&&
		( arg_beta_bridge_1.type != beta_bridge_type::ANTI_PARALLEL || arg_beta_bridge_1.partner_idx > arg_beta_bridge_2.partner_idx )
		&&
		( difference( arg_beta_bridge_1.partner_idx, arg_beta_bridge_2.partner_idx ) == 1 )
	);
}

/// \brief Set the contexts for the specified beta-bridges
void cath::sec::detail::set_bridges_contexts(beta_bridge_vec_vec &arg_beta_bridges ///< The beta-bridges to update
                                             ) {
	for (const size_t &beta_bridges_ctr : irange( 1_z, max( 2_z, arg_beta_bridges.size() ) - 1_z ) ) {
		const beta_bridge_vec &prev_res_bridges = arg_beta_bridges[ beta_bridges_ctr - 1 ];
		beta_bridge_vec       &this_res_bridges = arg_beta_bridges[ beta_bridges_ctr     ];
		const beta_bridge_vec &next_res_bridges = arg_beta_bridges[ beta_bridges_ctr + 1 ];

		for (beta_bridge &this_res_bridge : this_res_bridges) {
			const bool has_neighbour_bridge = (
				any_of( prev_res_bridges, [&] (const beta_bridge &prev_res_bridge) {
					return are_consecutive_bridges_in_sheet( prev_res_bridge, this_res_bridge );
				} )
				||
				any_of( next_res_bridges, [&] (const beta_bridge &next_res_bridge) {
					return are_consecutive_bridges_in_sheet( this_res_bridge, next_res_bridge );
				} )
			);

			this_res_bridge.context = has_neighbour_bridge ? beta_bridge_context::IN_SHEET
			                                               : beta_bridge_context::LONE_BRIDGE;
		}
	}
}

/// \brief Set the contexts for a copy of the specified beta-bridges and return the copy
beta_bridge_vec_vec cath::sec::detail::set_bridges_contexts_copy(beta_bridge_vec_vec arg_beta_bridges ///< The original beta-bridges
                                                                 ) {
	set_bridges_contexts( arg_beta_bridges );
	return arg_beta_bridges;
}

/// \brief Calculate the sec_struc_type values for the specified bifur_hbond_list
///
/// \relates bifur_hbond_list
sec_struc_type_vec cath::sec::calc_sec_strucs(const bifur_hbond_list &arg_bifur_hbond_list_raw ///< The bifur_hbond_list to query
                                              ) {
	const auto arg_bifur_hbond_list = remove_not_bondy_enough_copy( arg_bifur_hbond_list_raw );

	const auto beta_bridges = set_bridges_contexts_copy( transform_build<beta_bridge_vec_vec>(
		irange( 0_z, arg_bifur_hbond_list.size() ),
		[&] (const size_t &x) {
			return has_beta_bridge( arg_bifur_hbond_list, x );
		}
	) );

	// for (const size_t &beta_bridges_ctr : irange( 0_z, beta_bridges.size() ) ) {
	// 	std::cerr << beta_bridges_ctr << "\t" << join ( beta_bridges[ beta_bridges_ctr ] | lexical_casted<string>(), ", " ) << "\n";
	// }

	sec_struc_type_vec results( arg_bifur_hbond_list.size(), sec_struc_type::COIL );
	for (const size_t &beta_bridge_idx_1 : irange( 0_z, beta_bridges.size() ) ) {
		for (const size_t &beta_bridge_idx_2 : irange( beta_bridge_idx_1 + 1, min( beta_bridges.size(), beta_bridge_idx_1 + 3 ) ) ) {
			for (const auto &bridge_1 : beta_bridges[ beta_bridge_idx_1 ] ) {
				for (const auto &bridge_2 : beta_bridges[ beta_bridge_idx_2 ] ) {
					if ( is_beta_bulge( bridge_1, beta_bridge_idx_1, bridge_2, beta_bridge_idx_2 ) ) {
						// std::cerr << "Adding beta bulge from " << beta_bridge_idx_1 << " to " << beta_bridge_idx_2 << "\n";
						add_beta_bulge( results, beta_bridge_idx_1, beta_bridge_idx_2 );
						// std::cerr << "Adding beta bulge from " << bridge_1.partner_idx << " to " << bridge_2.partner_idx << "\n";
						add_beta_bulge( results, bridge_1.partner_idx, bridge_2.partner_idx );
					}
				}
			}
		}
	}

	// Label the residues
	for (const size_t &bifur_bond_ctr : irange( 0_z, arg_bifur_hbond_list.size() ) ) {
		// std::cerr << bifur_bond_ctr << "\t" << arg_bifur_hbond_list [ bifur_bond_ctr ] << "\n";

		// If this residue is part of a 4-helix (and not a 5-helix), label with alpha-helix
		if ( is_n_helix( arg_bifur_hbond_list, bifur_bond_ctr, 4 ) && ! is_n_helix( arg_bifur_hbond_list, bifur_bond_ctr, 5 ) ) {
			results[ bifur_bond_ctr ] = sec_struc_type::ALPHA_HELIX;
		}
		// If this residue has any beta-bridges that are in beta-sheets, label with beta-strand
		if ( any_of( beta_bridges[ bifur_bond_ctr ], [] (const beta_bridge &x) { return x.context == beta_bridge_context::IN_SHEET; } ) ) {
			results[ bifur_bond_ctr ] = sec_struc_type::BETA_STRAND;
		}
	}

	return results;
}

/// \brief Calculate the sec_struc_type values for the specified pdb
///
/// \relates pdb
sec_struc_type_vec cath::sec::calc_sec_strucs_of_pdb__recalc_backbone_residues(const pdb             &arg_pdb,   ///< The pdb to query
                                                                               const ostream_ref_opt &arg_stderr ///< An optional reference to an ostream to which any logging should be performed
                                                                               ) {
	return calc_sec_strucs( dssp_hbond_calc::calc_bifur_hbonds_of_pdb__recalc_backbone_residues(
		arg_pdb,
		arg_stderr
	) );
}

/// \brief Calculate the sec_struc_type values for the specified pdb
///
/// \relates pdb
sec_struc_type_vec cath::sec::calc_sec_strucs_of_backbone_complete_pdb(const pdb &arg_pdb ///< The pdb to query
                                                                       ) {
	return calc_sec_strucs( dssp_hbond_calc::calc_bifur_hbonds_of_backbone_complete_pdb(
		arg_pdb
	) );
}

/// \brief Calculate the sec_struc_type values for the specified protein
///
/// \relates protein
sec_struc_type_vec cath::sec::get_sec_strucs(const protein &arg_protein ///< The protein to query
                                             ) {
	return transform_build<sec_struc_type_vec>(
		arg_protein,
		[] (const residue &x) {
			return x.get_sec_struc_type();
		}
	);
}
