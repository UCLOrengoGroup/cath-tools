/// \file
/// \brief The data_dirs_options_block class definitions

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

#include "strucs_context.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/count_if.hpp>

#include "chopping/region/region.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/debug_numeric_cast.hpp"
#include "file/pdb/pdb.hpp"

using namespace cath::chop;
using namespace cath::common;
using namespace cath::file;

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::range::count_if;
using std::string;

/// \brief Restrict the PDBs in the specified strucs_context according to its regions
///
/// \relates strucs_context
void cath::file::restrict_pdbs(strucs_context &arg_strucs_context ///< The strucs_context to modify
                               ) {
	arg_strucs_context.set_pdbs( get_restricted_pdbs( arg_strucs_context ) );
}

/// \brief Make a copy of the specified strucs_context in which the PDBs are restricted according to its regions
///
/// \relates strucs_context
strucs_context cath::file::restrict_pdbs_copy(strucs_context arg_strucs_context ///< The source strucs_context
                                              ) {
	restrict_pdbs( arg_strucs_context );
	return arg_strucs_context;
}

/// \brief Get the number of entries in the specified strucs_context with regions specified
///
/// \relates strucs_context
size_t cath::file::get_num_regions_set(const strucs_context &arg_strucs_context ///< The strucs_context to query
                                       ) {
	return debug_numeric_cast<size_t>( count_if(
		arg_strucs_context.get_regions(),
		[] (const region_vec_opt &x) { return static_cast<bool>( x ); }
	) );
}

/// \brief Generate a string describing the specified strucs_context
///
/// \relates strucs_context
string cath::file::to_string(const strucs_context &arg_strucs_context ///< The strucs_context to describe
                             ) {
	return "strucs_context[" + join(
		indices( get_num_entries( arg_strucs_context ) )
			| transformed( [&] (const size_t &struc_context_idx) {
				using std::to_string;

				const pdb            &the_pdb         = arg_strucs_context.get_pdbs      () [ struc_context_idx ];
				const name_set       &the_name_set    = arg_strucs_context.get_name_sets () [ struc_context_idx ];
				const region_vec_opt &the_regions_opt = arg_strucs_context.get_regions   () [ struc_context_idx ];
				return "pdb["
					+ to_string( the_pdb.get_num_residues() )
					+ " residues]+"
					+ to_string( the_name_set )
					+ (
						the_regions_opt
						? "+" + join(
							*the_regions_opt
								| transformed( [] (const region &x) {
									return to_string( x );
								} ),
							","
						)
						: ""
					);
			} ),
		"; "
	) + "]";
}

/// \brief Get a copy of the PDBs in the specified strucs_context, restricted to its regions
///
/// \relates strucs_context
pdb_list cath::file::get_restricted_pdbs(const strucs_context &arg_strucs_context ///< The strucs_context to query
                                         ) {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return pdb_list{
		transform_build<pdb_vec>(
			arg_strucs_context.get_pdbs(),
			arg_strucs_context.get_regions(),
			[] (const pdb &the_pdb, const region_vec_opt &the_regions) {
				return get_regions_limited_pdb( the_regions, the_pdb );
			}
		)
	};
}
