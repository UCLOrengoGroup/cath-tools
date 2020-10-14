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
#include <boost/algorithm/cxx11/none_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/range/algorithm/binary_search.hpp>
#include <boost/range/algorithm/remove_if.hpp>
#include <boost/range/algorithm/upper_bound.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/append.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
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
using boost::algorithm::none_of;
using boost::format;
using boost::irange;
using boost::make_optional;
using boost::none;
using boost::optional;
using boost::range::binary_search;
using boost::range::upper_bound;
using boost::remove_if;
using std::max;
using std::min;
using std::ostream;
using std::string;

constexpr size_t sec_struc_consts::MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE;
constexpr size_t sec_struc_consts::BETA_BULGE_MAX_DIFF_SOURCE;
constexpr size_t sec_struc_consts::BETA_BULGE_MAX_DIFF_DEST;
constexpr size_t sec_struc_consts::DEFAULT_HELIX_N;
constexpr beta_bridge_context beta_bridge::DEFAULT_CONTEXT;

/// \brief Generate a string describing the specified beta_bridge_type
///
/// \relates beta_bridge_type
string cath::sec::detail::to_string(const beta_bridge_type &prm_beta_bridge_type ///< The beta_bridge_type to describe
                                    ) {
	switch ( prm_beta_bridge_type ) {
		case ( beta_bridge_type::PARALLEL      ) : { return "PARALLEL"; };
		case ( beta_bridge_type::ANTI_PARALLEL ) : { return "ANTIPARA"; };
	}
	BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of beta_bridge_type not recognised whilst converting to_string()"));
}

/// \brief Insert a description of the specified beta_bridge_type into the specified ostream
///
/// \relates beta_bridge_type
ostream & cath::sec::detail::operator<<(ostream                &prm_os,              ///< The ostream into which the description should be inserted
                                        const beta_bridge_type &prm_beta_bridge_type ///< The beta_bridge_type to describe
                                        ) {
	prm_os << to_string( prm_beta_bridge_type );
	return prm_os;
}

/// \brief Generate a string describing the specified beta_bridge_context
///
/// \relates beta_bridge_context
string cath::sec::detail::to_string(const beta_bridge_context &prm_beta_bridge_context ///< The beta_bridge_context to describe
                                    ) {
	switch ( prm_beta_bridge_context ) {
		case ( beta_bridge_context::LONE_BRIDGE ) : { return "BRDGE"; };
		case ( beta_bridge_context::IN_SHEET    ) : { return "SHEET"; };
	}
	BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of beta_bridge_context not recognised whilst converting to_string()"));
}

/// \brief Insert a description of the specified beta_bridge_context into the specified ostream
///
/// \relates beta_bridge_context
ostream & cath::sec::detail::operator<<(ostream                   &prm_os,                 ///< The ostream into which the description should be inserted
                                        const beta_bridge_context &prm_beta_bridge_context ///< The beta_bridge_context to describe
                                        ) {
	prm_os << to_string( prm_beta_bridge_context );
	return prm_os;
}

/// \brief Return whether the two specified beta_bridges are identical
///
/// \relates beta_bridge
inline bool cath::sec::detail::operator==(const beta_bridge &prm_lhs, ///< The first  beta_bridge to compare
                                          const beta_bridge &prm_rhs  ///< The second beta_bridge to compare
                                          ) {
	return (
		prm_lhs.partner_idx == prm_rhs.partner_idx
		&&
		prm_lhs.type        == prm_rhs.type
	);
}

/// \brief Generate a string describing the specified beta_bridge
///
/// \relates beta_bridge
string cath::sec::detail::to_string(const beta_bridge &prm_beta_bridge ///< The beta_bridge to describe
                                    ) {
	return "beta_bridge[partner:"
		+ ( format( "%5g") % prm_beta_bridge.partner_idx ).str()
		+ ", type:"
		+ to_string( prm_beta_bridge.type )
		+ ", context:"
		+ to_string( prm_beta_bridge.context )
		+ "]";
}

/// \brief Insert a description of the specified beta_bridge into the specified ostream
///
/// \relates beta_bridge
ostream & cath::sec::detail::operator<<(ostream           &prm_os,         ///< The ostream into which the description should be inserted
                                        const beta_bridge &prm_beta_bridge ///< The beta_bridge to describe
                                        ) {
	prm_os << to_string( prm_beta_bridge );
	return prm_os;
}

/// \brief Return whether the specified helix_category value means that the residue in question is helix-bonded to an earlier residue
///
/// \relates helix_category
bool cath::sec::detail::is_bonded_to_earlier(const helix_category &prm_helix_cat ///<  The helix_category to query
                                             ) {
	return ( prm_helix_cat == helix_category::BONDED_TO_BOTH || prm_helix_cat == helix_category::BONDED_TO_EARLIER_ONLY );
}

/// \brief Return whether the specified helix_category value means that the residue in question is helix-bonded to a later residue
///
/// \relates helix_category
bool cath::sec::detail::is_bonded_to_later(const helix_category &prm_helix_cat ///< The helix_category to query
                                           ) {
	return ( prm_helix_cat == helix_category::BONDED_TO_BOTH || prm_helix_cat == helix_category::BONDED_TO_LATER_ONLY );
}


/// \brief Return whether the specified hbond_half_opt_pair is bonded to the residue at the specified index
bool cath::sec::detail::is_bonded_to(const hbond_half_opt_pair &prm_bound_pair, ///< The hbond_half_opt_pair to examine
                                     const size_t              &prm_dest_index  ///< The index of the residue to which we're interested in potential bonds
                                     ) {
	return (
		is_bondy_enough( prm_bound_pair.first )
		&&
		(
			( prm_bound_pair.first->index == prm_dest_index )
			||
			(
				is_bondy_enough( prm_bound_pair.second )
				&&
				prm_bound_pair.second->index == prm_dest_index
			)
		)
	);
}

/// \brief Return whether there is an adequate bond between the NH of the residue at the first specified
///        index and the CO of the residue at the second, and in that direction.
bool cath::sec::detail::are_strictly_nh_to_co_bonded(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                     const size_t           &prm_nh_res_idx,       ///< The index of the residue at the NH side of the required bond
                                                     const size_t           &prm_co_res_idx        ///< The index of the residue at the CO side of the required bond
                                                     ) {
	return ( is_bonded_to( prm_bifur_hbond_list[ prm_nh_res_idx ].get_bound_pair_for_this_nh(), prm_co_res_idx ) );
}

/// \brief Return whether there is an adequate bond between the NH of the residue at the first specified
///        index and the CO of the residue at the second
///
/// This will accept the bond being found at either side alone
bool cath::sec::detail::are_nh_to_co_bonded(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                            const size_t           &prm_nh_res_idx,       ///< The index of the residue at the NH side of the required bond
                                            const size_t           &prm_co_res_idx        ///< The index of the residue at the CO side of the required bond
                                            ) {
	return (
		is_bonded_to( prm_bifur_hbond_list[ prm_nh_res_idx ].get_bound_pair_for_this_nh(), prm_co_res_idx )
		||
		is_bonded_to( prm_bifur_hbond_list[ prm_co_res_idx ].get_bound_pair_for_this_co(), prm_nh_res_idx )
	);
}

/// \brief Return whether there is an adequate bond between the CO of the residue at the first specified
///        index and the NH of the residue at the second
///
/// This will accept the bond being found at either side alone
bool cath::sec::detail::are_co_to_nh_bonded(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                            const size_t           &prm_co_res_idx,       ///< The index of the residue at the CO side of the required bond
                                            const size_t           &prm_nh_res_idx        ///< The index of the residue at the NH side of the required bond
                                            ) {
	return (
		is_bonded_to( prm_bifur_hbond_list[ prm_co_res_idx ].get_bound_pair_for_this_co(), prm_nh_res_idx )
		||
		is_bonded_to( prm_bifur_hbond_list[ prm_nh_res_idx ].get_bound_pair_for_this_nh(), prm_co_res_idx )
	);
}

/// \brief Get the n-helix category of the specified bifur_hbond at the specified index
///
/// \pre `prm_break_indices` must be sorted in ascending order
optional<helix_category> cath::sec::detail::n_helix_cat(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond to query
                                                        const size_vec         &prm_break_indices,    ///< A list of the residues that are preceded by a chain break in ascending order
                                                        const size_t           &prm_index,            ///< The index of the bifur_hbond
                                                        const size_t           &prm_helix_num         ///< "n", the difference between bonded residues
                                                        ) {
	const bool bonded_to_later   = (
		prm_index + prm_helix_num < prm_bifur_hbond_list.size()
		&&
		are_strictly_nh_to_co_bonded( prm_bifur_hbond_list, prm_index + prm_helix_num, prm_index                 )
		&&
		! indices_straddle_break    ( prm_break_indices,    prm_index + prm_helix_num, prm_index                 )
	);
	const bool bonded_to_earlier = (
		prm_index >= prm_helix_num
		&&
		are_strictly_nh_to_co_bonded( prm_bifur_hbond_list, prm_index,                 prm_index - prm_helix_num )
		&&
		! indices_straddle_break    ( prm_break_indices,    prm_index,                 prm_index - prm_helix_num )
	);

	if ( bonded_to_later ) {
		return { bonded_to_earlier ? helix_category::BONDED_TO_BOTH : helix_category::BONDED_TO_LATER_ONLY };
	}
	return make_optional( bonded_to_earlier, helix_category::BONDED_TO_EARLIER_ONLY );
}

/// \brief Return whether the residue at the specified is n-helix bonded to the relevant later residue
///
/// \pre `prm_break_indices` must be sorted in ascending order
bool cath::sec::detail::is_n_helix_bonded_to_later(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond to query
                                                   const size_vec         &prm_break_indices,    ///< A list of the residues that are preceded by a chain break in ascending order
                                                   const size_t           &prm_index,            ///< The index of the residue to query
                                                   const size_t           &prm_helix_num         ///< "n", the difference between bonded residues
                                                   ) {
	const auto cat_opt = n_helix_cat( prm_bifur_hbond_list, prm_break_indices, prm_index, prm_helix_num );
	return ( cat_opt && is_bonded_to_later( *cat_opt ) );
}

/// \brief Return whether the residue at the specified index could potentially start n-helix according to the specified bifur_hbond_list
///        (not considering clashes with other types of helix)
///
/// \pre `prm_break_indices` must be sorted in ascending order
bool cath::sec::detail::could_start_n_helix(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                            const size_vec         &prm_break_indices,    ///< A list of the residues that are preceded by a chain break in ascending order
                                            const size_t           &prm_index,            ///< The index of the residue in question
                                            const size_t           &prm_helix_num         ///< "n", the difference between bonded residues
                                            ) {
	if ( prm_index == 0 ) {
		return false;
	}
	return (
		is_n_helix_bonded_to_later( prm_bifur_hbond_list, prm_break_indices, prm_index - 1, prm_helix_num )
		&&
		is_n_helix_bonded_to_later( prm_bifur_hbond_list, prm_break_indices, prm_index,     prm_helix_num )
	);
}

/// \brief Return whether the residue at the specified index is at the start of a 3-helix according to the specified bifur_hbond_list
///
/// \pre `prm_break_indices` must be sorted in ascending order
bool cath::sec::detail::starts_3_helix(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                       const size_vec         &prm_break_indices,    ///< A list of the residues that are preceded by a chain break in ascending order
                                       const size_t           &prm_index             ///< The index of the residue in question
                                       ) {
	return (
		! prm_bifur_hbond_list.empty()
		&&
		could_start_n_helix( prm_bifur_hbond_list, prm_break_indices, prm_index, 3 )
		&&
		none_of(
			irange( max( 3_z, prm_index ) - 3_z, min( prm_index + 3_z, prm_bifur_hbond_list.size() ) ),
			[&] (const size_t &x) { return could_start_n_helix( prm_bifur_hbond_list, prm_break_indices, x, 4 ); }
		)
	);

}

/// \brief Return whether the residue at the specified index is at the start of a 5-helix according to the specified bifur_hbond_list
///
/// \pre `prm_break_indices` must be sorted in ascending order
bool cath::sec::detail::starts_5_helix(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                       const size_vec         &prm_break_indices,    ///< A list of the residues that are preceded by a chain break in ascending order
                                       const size_t           &prm_index             ///< The index of the residue in question
                                       ) {
	return (
		! prm_bifur_hbond_list.empty()
		&&
		could_start_n_helix( prm_bifur_hbond_list, prm_break_indices, prm_index, 5 )
		&&
		none_of(
			irange( max( prm_index, 2_z) - 2_z, min( prm_index + 5_z, prm_bifur_hbond_list.size() ) ),
			[&] (const size_t &x) { return starts_3_helix( prm_bifur_hbond_list, prm_break_indices, x ); }
		)
	);

}

/// \brief Return whether the residue at the specified index is within a 5-helix according to the specified bifur_hbond_list
///
/// \pre `prm_break_indices` must be sorted in ascending order
bool cath::sec::detail::in_5_helix(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                   const size_vec         &prm_break_indices,    ///< A list of the residues that are preceded by a chain break in ascending order
                                   const size_t           &prm_index             ///< The index of the residue in question
                                   ) {
	return any_of(
		irange( max( 5_z, prm_index ) - 4_z, prm_index + 1_z ) | reversed,
		[&] (const size_t &x) { return starts_5_helix( prm_bifur_hbond_list, prm_break_indices, x ); }
	);
}

/// \brief Return whether the residue at the specified index is within a four-helix but not five-helix according to the specified bifur_hbond_list
///
/// \pre `prm_break_indices` must be sorted in ascending order
bool cath::sec::detail::is_in_4_helix_not_conflicting_with_5_helix(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                   const size_vec         &prm_break_indices,    ///< A list of the residues that are preceded by a chain break in ascending order
                                                                   const size_t           &prm_index             ///< The index of the residue in question
                                                                   ) {
	return (
		any_of(
			irange( max( 4_z, prm_index ) - 3_z, prm_index + 1_z ) | reversed,
			[&] (const size_t &x) { return could_start_n_helix( prm_bifur_hbond_list, prm_break_indices, x, 4_z ); }
		)
		&&
		! in_5_helix( prm_bifur_hbond_list, prm_break_indices, prm_index )
	);
}

/// \brief Return whether the residue at the specified index is within a four-helix according to the specified bifur_hbond_list
///
/// \pre `prm_break_indices` must be sorted in ascending order
bool cath::sec::detail::is_in_n_helix(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                      const size_vec         &prm_break_indices,    ///< A list of the residues that are preceded by a chain break in ascending order
                                      const size_t           &prm_index,            ///< The index of the residue in question
                                      const size_t           &prm_helix_num         ///< "n", the difference between bonded residues
                                      ) {
	return any_of(
		irange( max( prm_helix_num, prm_index ) - ( prm_helix_num - 1_z ), prm_index + 1_z ) | reversed,
		[&] (const size_t &x) { return could_start_n_helix( prm_bifur_hbond_list, prm_break_indices, x, prm_helix_num ); }
	);
}

/// \brief Return whether the specified index is within range to be the middle of a beta triplet
bool cath::sec::detail::beta_index_in_range(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list providing the context
                                            const size_t           &prm_index             ///< The index to check
                                            ) {
	return (
		( prm_index > 0 )
		&&
		( prm_index + 1 < prm_bifur_hbond_list.size() )
	);
}

/// \brief Return whether there are parallel-beta-bridge bonds between the specified indices
///        that bond to the source residue itself
beta_bridge_opt cath::sec::detail::has_parallel_beta_bridge_bonds_to_src(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                         const size_t           &prm_src_index,        ///< The index of the source residue to query
                                                                         const size_t           &prm_dest_index        ///< The index of the destination residue to query
                                                                         ) {
	return make_optional(
		(
			difference( prm_src_index, prm_dest_index ) >= sec_struc_consts::MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE
			&&
			beta_index_in_range( prm_bifur_hbond_list, prm_src_index  )
			&&
			beta_index_in_range( prm_bifur_hbond_list, prm_dest_index )
			&&
			are_strictly_nh_to_co_bonded( prm_bifur_hbond_list, prm_src_index,      prm_dest_index - 1 )
			&&
			are_strictly_nh_to_co_bonded( prm_bifur_hbond_list, prm_dest_index + 1, prm_src_index      )
		),
		beta_bridge{ prm_dest_index, beta_bridge_type::PARALLEL }
	);
}

/// \brief Return whether there are parallel-beta-bridge bonds between the specified indices
///        that bond to the residues straddling the source residue
beta_bridge_opt cath::sec::detail::has_parallel_beta_bridge_bonds_straddling_src(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                                 const size_t           &prm_src_index,        ///< The index of the source residue to query
                                                                                 const size_t           &prm_dest_index        ///< The index of the destination residue to query
                                                                                 ) {
	return make_optional(
		(
			difference( prm_src_index, prm_dest_index ) >= sec_struc_consts::MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE
			&&
			beta_index_in_range( prm_bifur_hbond_list, prm_src_index  )
			&&
			beta_index_in_range( prm_bifur_hbond_list, prm_dest_index )
			&&
			are_strictly_nh_to_co_bonded( prm_bifur_hbond_list, prm_dest_index,    prm_src_index - 1 )
			&&
			are_strictly_nh_to_co_bonded( prm_bifur_hbond_list, prm_src_index + 1, prm_dest_index    )
		),
		beta_bridge{ prm_dest_index, beta_bridge_type::PARALLEL }
	);
}

/// \brief Return whether there are antiparallel-beta-bridge bonds between the specified indices
///        that bond to the source residue itself
beta_bridge_opt cath::sec::detail::has_antiparallel_beta_bridge_bonds_to_src(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                             const size_t           &prm_src_index,        ///< The index of the source residue to query
                                                                             const size_t           &prm_dest_index        ///< The index of the destination residue to query
                                                                             ) {
	return make_optional(
		(
			difference( prm_src_index, prm_dest_index ) >= sec_struc_consts::MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE
			&&
			prm_src_index  + 1 < prm_bifur_hbond_list.size()
			&&
			prm_dest_index + 1 < prm_bifur_hbond_list.size()
			&&
			are_strictly_nh_to_co_bonded( prm_bifur_hbond_list, prm_src_index,  prm_dest_index )
			&&
			are_strictly_nh_to_co_bonded( prm_bifur_hbond_list, prm_dest_index, prm_src_index  )
		),
		beta_bridge{ prm_dest_index, beta_bridge_type::ANTI_PARALLEL }
	);
}

/// \brief Return whether there are antiparallel-beta-bridge bonds between the specified indices
///        that bond to the residues straddling the source residue
beta_bridge_opt cath::sec::detail::has_antiparallel_beta_bridge_bonds_straddling_src(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                                     const size_t           &prm_src_index,        ///< The index of the source residue to query
                                                                                     const size_t           &prm_dest_index        ///< The index of the destination residue to query
                                                                                     ) {
	return make_optional(
		(
			difference( prm_src_index, prm_dest_index ) >= sec_struc_consts::MIN_ALLOWABLE_RES_DIFF_FOR_BETA_BRIDGE
			&&
			beta_index_in_range( prm_bifur_hbond_list, prm_src_index  )
			&&
			beta_index_in_range( prm_bifur_hbond_list, prm_dest_index )
			&&
			are_strictly_nh_to_co_bonded( prm_bifur_hbond_list, prm_dest_index + 1, prm_src_index  - 1 )
			&&
			are_strictly_nh_to_co_bonded( prm_bifur_hbond_list, prm_src_index  + 1, prm_dest_index - 1 )
		),
		beta_bridge{ prm_dest_index, beta_bridge_type::ANTI_PARALLEL }
	);
}

/// \brief Return whether the residue at the specified index is part of a parallel beta bridge according to the specified bifur_hbond_list
beta_bridge_vec cath::sec::detail::has_parallel_beta_bridge(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                            const size_t           &prm_index             ///< The index to query
                                                            ) {
	if ( ! beta_index_in_range( prm_bifur_hbond_list, prm_index ) ) {
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
	const auto to_src_fn    = [&] (const size_t &x) { cnsdr_res_fn( has_parallel_beta_bridge_bonds_to_src        ( prm_bifur_hbond_list, prm_index, x ) ); };
	const auto strdl_src_fn = [&] (const size_t &x) { cnsdr_res_fn( has_parallel_beta_bridge_bonds_straddling_src( prm_bifur_hbond_list, prm_index, x ) ); };

	const auto prev_co_fn   = [&] (const auto   &x) { if ( is_bondy_enough( x )                     ) { strdl_src_fn( x->index     ); } };
	const auto this_fn      = [&] (const auto   &x) { if ( is_bondy_enough( x ) && ( x->index > 1 ) ) { to_src_fn   ( x->index - 1 ); } };
	const auto next_nh_fn   = [&] (const auto   &x) { if ( is_bondy_enough( x )                     ) { strdl_src_fn( x->index     ); } };

	const auto &prev_co = prm_bifur_hbond_list[ prm_index - 1 ].get_bound_pair_for_this_co();
	const auto &this_co = prm_bifur_hbond_list[ prm_index     ].get_bound_pair_for_this_co();
	const auto &next_nh = prm_bifur_hbond_list[ prm_index + 1 ].get_bound_pair_for_this_nh();

	prev_co_fn( prev_co.first  );
	prev_co_fn( prev_co.second );
	this_fn   ( this_co.first  );
	this_fn   ( this_co.second );
	next_nh_fn( next_nh.first  );
	next_nh_fn( next_nh.second );
	return results;
}

/// \brief Return whether the residue at the specified index is part of an antiparallel beta bridge according to the specified bifur_hbond_list
beta_bridge_vec cath::sec::detail::has_antiparallel_beta_bridge(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                                const size_t           &prm_index             ///< The index to query
                                                                ) {
	if ( ! beta_index_in_range( prm_bifur_hbond_list, prm_index ) ) {
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

	const auto sz = prm_bifur_hbond_list.size();

	const auto to_src_fn    = [&] (const size_t &x) { cnsdr_res_fn( has_antiparallel_beta_bridge_bonds_to_src        ( prm_bifur_hbond_list, prm_index, x ) ); };
	const auto strdl_src_fn = [&] (const size_t &x) { cnsdr_res_fn( has_antiparallel_beta_bridge_bonds_straddling_src( prm_bifur_hbond_list, prm_index, x ) ); };

	/// \todo Come GCCs all > 5.2.1, use a generic lambda here (causes ICE in GCC 5.2.1)
	const auto prev_co_fn   = [&] (const hbond_half_opt &x) { if ( is_bondy_enough( x ) && (   x->index       > 1  ) ) { strdl_src_fn( x->index - 1 ); } };
	const auto this_fn      = [&] (const hbond_half_opt &x) { if ( is_bondy_enough( x )                              ) { to_src_fn   ( x->index     ); } };
	const auto next_nh_fn   = [&] (const hbond_half_opt &x) { if ( is_bondy_enough( x ) && ( ( x->index + 1 ) < sz ) ) { strdl_src_fn( x->index + 1 ); } };

	const auto &prev_co = prm_bifur_hbond_list[ prm_index - 1 ].get_bound_pair_for_this_co();
	const auto &this_nh = prm_bifur_hbond_list[ prm_index     ].get_bound_pair_for_this_nh();
	const auto &this_co = prm_bifur_hbond_list[ prm_index     ].get_bound_pair_for_this_co();
	const auto &next_nh = prm_bifur_hbond_list[ prm_index + 1 ].get_bound_pair_for_this_nh();

	prev_co_fn( prev_co.first  );
	prev_co_fn( prev_co.second );
	this_fn   ( this_nh.first  );
	this_fn   ( this_nh.second );
	this_fn   ( this_co.first  );
	this_fn   ( this_co.second );
	next_nh_fn( next_nh.first  );
	next_nh_fn( next_nh.second );
	return results;
}

/// \brief Return whether the specified index is part of a beta bridge according to the specified bifur_hbond_list
beta_bridge_vec cath::sec::detail::has_beta_bridge(const bifur_hbond_list &prm_bifur_hbond_list, ///< The bifur_hbond_list to query
                                                   const size_t           &prm_index             ///< The index to query
                                                   ) {
	return append_copy(
		has_parallel_beta_bridge    ( prm_bifur_hbond_list, prm_index ),
		has_antiparallel_beta_bridge( prm_bifur_hbond_list, prm_index )
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
bool cath::sec::detail::is_beta_bulge(const beta_bridge &prm_bridge_1, ///< The first  beta-bridge
                                      const size_t      &prm_index_1,  ///< The source index of the first  beta-bridge
                                      const beta_bridge &prm_bridge_2, ///< The second beta-bridge
                                      const size_t      &prm_index_2   ///< The source index of the second beta-bridge
                                      ) {
	const size_t smaller_1       = min( prm_index_1, prm_bridge_1.partner_idx );
	const size_t smaller_2       = min( prm_index_2, prm_bridge_2.partner_idx );
	const size_t smaller_diff    = difference( smaller_1, smaller_2 );
	const size_t src_index_diff  = difference( prm_index_1,              prm_index_2              );
	const size_t dest_index_diff = difference( prm_bridge_1.partner_idx, prm_bridge_2.partner_idx );
	const bool   src_ascending   = ( prm_index_2              > prm_index_1              );
	const bool   dest_ascending  = ( prm_bridge_2.partner_idx > prm_bridge_1.partner_idx );
	const bool   dirns_match     = ( src_ascending == dest_ascending );
	return (
		prm_bridge_1.type == prm_bridge_2.type
		&&
		clamp( src_index_diff,  0, sec_struc_consts::BETA_BULGE_MAX_DIFF_SOURCE ) == src_index_diff
		&&
		clamp( dest_index_diff, 0, sec_struc_consts::BETA_BULGE_MAX_DIFF_DEST   ) == dest_index_diff
		&&
		smaller_diff > 0
		&&
		max( src_index_diff, dest_index_diff ) > 1
		&&
		(
			( prm_bridge_1.type == beta_bridge_type::PARALLEL      &&   dirns_match )
			||
			( prm_bridge_2.type == beta_bridge_type::ANTI_PARALLEL && ! dirns_match )
		)
	);
}

/// \relates Add a beta-bulge between the specified residues in the specified sec_struc_types
///
/// Note: the from/to residues needn't necessarily be in ascending order
void cath::sec::detail::add_beta_bulge(sec_struc_type_vec &prm_secs, ///< The vector of sec_struc_type values to update
                                       const size_t       &prm_from, ///< The index of the residue at one end of the bulge
                                       const size_t       &prm_to    ///< The index of the residue at the other end of the bulge
                                       ) {
	for (const size_t &x : irange( min( prm_from, prm_to ), max( prm_from, prm_to ) + 1_z ) ) {
		prm_secs[ x ] = sec_struc_type::BETA_STRAND;
	}
}

/// \brief Return whether the two specified beta bridges are consecutive bridges in a beta sheet
///        (assuming the first's residue immediately follows the second's)
bool cath::sec::detail::are_consecutive_bridges_in_sheet(const beta_bridge &prm_beta_bridge_1, ///< The first  bridge to query
                                                         const beta_bridge &prm_beta_bridge_2  ///< The second bridge to query
                                                         ) {
	return (
		( prm_beta_bridge_1.type == prm_beta_bridge_2.type )
		&&
		( prm_beta_bridge_1.type != beta_bridge_type::PARALLEL      || prm_beta_bridge_1.partner_idx < prm_beta_bridge_2.partner_idx )
		&&
		( prm_beta_bridge_1.type != beta_bridge_type::ANTI_PARALLEL || prm_beta_bridge_1.partner_idx > prm_beta_bridge_2.partner_idx )
		&&
		( difference( prm_beta_bridge_1.partner_idx, prm_beta_bridge_2.partner_idx ) == 1 )
	);
}

/// \brief Set the contexts for the specified beta-bridges
void cath::sec::detail::set_bridges_contexts(beta_bridge_vec_vec &prm_beta_bridges ///< The beta-bridges to update
                                             ) {
	for (const size_t &beta_bridges_ctr : irange( 1_z, max( 2_z, prm_beta_bridges.size() ) - 1_z ) ) {
		const beta_bridge_vec &prev_res_bridges = prm_beta_bridges[ beta_bridges_ctr - 1 ];
		beta_bridge_vec       &this_res_bridges = prm_beta_bridges[ beta_bridges_ctr     ];
		const beta_bridge_vec &next_res_bridges = prm_beta_bridges[ beta_bridges_ctr + 1 ];

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
beta_bridge_vec_vec cath::sec::detail::set_bridges_contexts_copy(beta_bridge_vec_vec prm_beta_bridges ///< The original beta-bridges
                                                                 ) {
	set_bridges_contexts( prm_beta_bridges );
	return prm_beta_bridges;
}

/// \brief Remove any of the specified bridges to/from any residues on either side of a chain break, according to the specified list
///
/// \pre `prm_break_indices` must be sorted in ascending order
void cath::sec::detail::remove_bridges_to_chain_break_residues(beta_bridge_vec_vec &prm_beta_bridges, ///< The original beta-bridges
                                                               const size_vec      &prm_break_indices ///< A list of the residues that are preceded by a chain break in ascending order
                                                               ) {
	// For each break index, remove all lists of beta bridges on either side of the chain break
	for (const size_t &x : prm_break_indices) {
		prm_beta_bridges[ x                   ].clear();
		prm_beta_bridges[ max( 1_z, x ) - 1_z ].clear();
	}

	// For each residue's worth of bridges...
	for (beta_bridge_vec &res_bridges : prm_beta_bridges) {
		// Remove any bridges to either side of any of the breaks
		res_bridges.erase(
			remove_if(
				res_bridges,
				[&] (const beta_bridge &x) {
					return (
						binary_search( prm_break_indices, x.partner_idx )
						||
						(
							x.partner_idx + 1 < prm_beta_bridges.size()
							&&
							binary_search( prm_break_indices, x.partner_idx + 1 )
						)
					);
				}
			),
			common::cend( res_bridges )
		);
	}
}

/// \brief Take a copy of the specified bridges, remove any to/from any residues on either side of a chain break, according to the specified list
///        and return the copy
///
/// \pre `prm_break_indices` must be sorted in ascending order
beta_bridge_vec_vec cath::sec::detail::remove_bridges_to_chain_break_residues_copy(beta_bridge_vec_vec  prm_beta_bridges, ///< The original beta-bridges
                                                                                   const size_vec      &prm_break_indices ///< A list of the residues that are preceded by a chain break in ascending order
                                                                                   ) {
	remove_bridges_to_chain_break_residues( prm_beta_bridges, prm_break_indices );
	return prm_beta_bridges;
}

/// \brief Return whether the specified indices straddle any of the breaks represented by the specified list
///
/// \pre `prm_break_indices` must be sorted in ascending order
bool cath::sec::detail::indices_straddle_break(const size_vec &prm_break_indices, ///< A list of the indices of the residues that are preceded by a chain break in ascending order
                                               const size_t   &prm_index_a,       ///< The first  index to examine
                                               const size_t   &prm_index_b        ///< The second index to examine
                                               ) {
	return (
		upper_bound( prm_break_indices, prm_index_a )
		!=
		upper_bound( prm_break_indices, prm_index_b )
	);
}

/// \brief Calculate the sec_struc_type values for the specified bifur_hbond_list
///
/// \pre `prm_break_indices` must be sorted in ascending order
///
/// \relates bifur_hbond_list
sec_struc_type_vec cath::sec::calc_sec_strucs(const bifur_hbond_list &prm_bifur_hbond_list_raw, ///< The bifur_hbond_list to query
                                              const size_vec         &prm_break_indices         ///< A list of the residues that are preceded by a chain break in ascending order
                                              ) {
	const auto prm_bifur_hbond_list = remove_not_bondy_enough_copy( prm_bifur_hbond_list_raw );

	const auto beta_bridges = set_bridges_contexts_copy( remove_bridges_to_chain_break_residues_copy(
		transform_build<beta_bridge_vec_vec>(
			indices( prm_bifur_hbond_list.size() ),
			[&] (const size_t &x) {
				return has_beta_bridge( prm_bifur_hbond_list, x );
			}
		),
		prm_break_indices
	) );

	// for (const size_t &beta_bridges_ctr : indices( beta_bridges.size() ) ) {
	// 	std::cerr << beta_bridges_ctr << "\t" << join ( beta_bridges[ beta_bridges_ctr ] | lexical_casted<string>(), ", " ) << "\n";
	// }

	sec_struc_type_vec results( prm_bifur_hbond_list.size(), sec_struc_type::COIL );
	for (const size_t &beta_bridge_idx_1 : indices( beta_bridges.size() ) ) {
		for (const size_t &beta_bridge_idx_2 : irange( beta_bridge_idx_1, min( beta_bridges.size(), beta_bridge_idx_1 + 3 ) ) ) {
			for (const auto &bridge_1 : beta_bridges[ beta_bridge_idx_1 ] ) {
				for (const auto &bridge_2 : beta_bridges[ beta_bridge_idx_2 ] ) {
					if ( beta_bridge_idx_1 == beta_bridge_idx_2 && bridge_1 == bridge_2 ) {
						continue;
					}
					const bool is_beta_bulge_decn = (
						! indices_straddle_break( prm_break_indices, beta_bridge_idx_1,    beta_bridge_idx_2    )
						&&
						! indices_straddle_break( prm_break_indices, bridge_1.partner_idx, bridge_2.partner_idx )
						&&
						is_beta_bulge( bridge_1, beta_bridge_idx_1, bridge_2, beta_bridge_idx_2 )
					);
					if ( is_beta_bulge_decn ) {
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
	for (const size_t &bifur_bond_ctr : indices( prm_bifur_hbond_list.size() ) ) {
		// std::cerr << bifur_bond_ctr << "\t" << prm_bifur_hbond_list [ bifur_bond_ctr ];
		// std::cerr <<  "  [";
		// std::cerr << ( is_n_helix_bonded_to_later( prm_bifur_hbond_list, prm_break_indices, bifur_bond_ctr, 4 ) ? "4"s : " "s );
		// std::cerr <<  "] [";
		// std::cerr << ( is_n_helix_bonded_to_later( prm_bifur_hbond_list, prm_break_indices, bifur_bond_ctr, 3 ) ? "3"s : " "s );
		// std::cerr <<  "] [";
		// std::cerr << ( is_n_helix_bonded_to_later( prm_bifur_hbond_list, prm_break_indices, bifur_bond_ctr, 5 ) ? "5"s : " "s );
		// std::cerr <<  "]\n";

		// If this residue is part of a 4-helix (and not a 5-helix), label with alpha-helix
		if ( is_in_4_helix_not_conflicting_with_5_helix( prm_bifur_hbond_list, prm_break_indices, bifur_bond_ctr ) ) {
			results[ bifur_bond_ctr ] = sec_struc_type::ALPHA_HELIX;
		}
		// If this residue has any beta-bridges that are in beta-sheets, label with beta-strand
		else if ( any_of( beta_bridges[ bifur_bond_ctr ], [] (const beta_bridge &x) { return x.context == beta_bridge_context::IN_SHEET; } ) ) {
			results[ bifur_bond_ctr ] = sec_struc_type::BETA_STRAND;
		}
	}

	return results;
}

/// \brief Calculate the sec_struc_type values for the specified pdb
///
/// \relates pdb
sec_struc_type_vec cath::sec::calc_sec_strucs_of_pdb__recalc_backbone_residues(const pdb             &prm_pdb,   ///< The pdb to query
                                                                               const ostream_ref_opt &prm_stderr ///< An optional reference to an ostream to which any logging should be performed
                                                                               ) {
	const auto backbone_pdb = backbone_complete_subset_of_pdb(
		prm_pdb,
		prm_stderr,
		dssp_skip_res_skipping::SKIP
	).first;
	return calc_sec_strucs(
		dssp_hbond_calc::calc_bifur_hbonds_of_backbone_complete_pdb( backbone_pdb ),
		indices_of_residues_following_chain_breaks( backbone_pdb )
	);
}

/// \brief Calculate the sec_struc_type values for the specified pdb
///
/// \relates pdb
sec_struc_type_vec cath::sec::calc_sec_strucs_of_backbone_complete_pdb(const pdb &prm_pdb ///< The pdb to query
                                                                       ) {
	return calc_sec_strucs(
		dssp_hbond_calc::calc_bifur_hbonds_of_backbone_complete_pdb( prm_pdb ),
		indices_of_residues_following_chain_breaks( prm_pdb )
	);
}

/// \brief Calculate the sec_struc_type values for the specified protein
///
/// \relates protein
sec_struc_type_vec cath::sec::get_sec_strucs(const protein &prm_protein ///< The protein to query
                                             ) {
	return transform_build<sec_struc_type_vec>(
		prm_protein,
		[] (const residue &x) {
			return x.get_sec_struc_type();
		}
	);
}
