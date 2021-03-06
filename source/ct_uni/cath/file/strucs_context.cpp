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

#include "cath/chopping/domain/domain.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/hash/hash_value_combine.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/strucs_context.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_list.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;

using ::boost::adaptors::transformed;
using ::boost::algorithm::join;
using ::boost::range::count_if;
using ::std::hash;
using ::std::nullopt;
using ::std::string;

/// \brief Build a strucs_context from the specified one that contains the backbone-complete subsets of the PDBs
///
/// \relates strucs_context
strucs_context cath::file::strucs_context_of_backbone_complete_subset_pdbs(const strucs_context  &prm_strucs_context, ///< The strucs_context to use as a source
                                                                           const ostream_ref_opt &prm_ostream         ///< An optional reference to an ostream to which any logging should be sent
                                                                           ) {
	return {
		pdb_list_of_backbone_complete_subset_pdbs(
			prm_strucs_context.get_pdbs(),
			prm_ostream
		),
		prm_strucs_context.get_name_sets(),
		prm_strucs_context.get_regions()
	};
}

/// \brief Build a strucs_context from the specified one that contains the backbone-complete subsets of the PDBs
///        and is restricted to the regions specified in the strucs_context
///
/// \relates strucs_context
strucs_context cath::file::strucs_context_of_backbone_complete_region_limited_subset_pdbs(const strucs_context     &prm_strucs_context, ///< The strucs_context to use as a source
                                                                                          const ostream_ref_opt    &prm_ostream         ///< An optional reference to an ostream to which any logging should be sent
                                                                                          ) {
	return {
		pdb_list_of_backbone_complete_region_limited_subset_pdbs(
			prm_strucs_context.get_pdbs(),
			prm_strucs_context.get_regions(),
			prm_ostream
		),
		prm_strucs_context.get_name_sets(),
		prm_strucs_context.get_regions()
	};
}

/// \brief Restrict the PDBs in the specified strucs_context according to its regions
///
/// \relates strucs_context
void cath::file::restrict_pdbs(strucs_context &prm_strucs_context ///< The strucs_context to modify
                               ) {
	prm_strucs_context.set_pdbs( get_restricted_pdbs( prm_strucs_context ) );
}

/// \brief Make a copy of the specified strucs_context in which the PDBs are restricted according to its regions
///
/// \relates strucs_context
strucs_context cath::file::restrict_pdbs_copy(strucs_context prm_strucs_context ///< The source strucs_context
                                              ) {
	restrict_pdbs( prm_strucs_context );
	return prm_strucs_context;
}

/// \brief Get the number of entries in the specified strucs_context with regions specified
///
/// \relates strucs_context
size_t cath::file::get_num_regions_set(const strucs_context &prm_strucs_context ///< The strucs_context to query
                                       ) {
	return debug_numeric_cast<size_t>( count_if(
		prm_strucs_context.get_regions(),
		[] (const region_vec_opt &x) { return static_cast<bool>( x ); }
	) );
}

/// \brief Generate a string describing the specified strucs_context
///
/// \relates strucs_context
string cath::file::to_string(const strucs_context &prm_strucs_context ///< The strucs_context to describe
                             ) {
	return "strucs_context[" + join(
		indices( size( prm_strucs_context ) )
			| transformed( [&] (const size_t &struc_context_idx) {
				using ::std::to_string;

				const pdb            &the_pdb         = prm_strucs_context.get_pdbs      () [ struc_context_idx ];
				const name_set       &the_name_set    = prm_strucs_context.get_name_sets () [ struc_context_idx ];
				const region_vec_opt &the_regions_opt = prm_strucs_context.get_regions   () [ struc_context_idx ];
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
pdb_list cath::file::get_restricted_pdbs(const strucs_context &prm_strucs_context ///< The strucs_context to query
                                         ) {
	return pdb_list{
		transform_build<pdb_vec>(
			prm_strucs_context.get_pdbs(),
			prm_strucs_context.get_regions(),
			[] (const pdb &the_pdb, const region_vec_opt &the_regions) {
				return get_regions_limited_pdb( the_regions, the_pdb );
			}
		)
	};
}

/// \brief Make a protein list from the PDBs and names in the specified strucs_context
///
/// \relates strucs_context
protein_list cath::file::build_protein_list(const strucs_context &prm_strucs_context ///< The strucs_context from which to build the protein_list
                                            ) {
	return build_protein_list_of_pdb_list_and_names(
		prm_strucs_context.get_pdbs(),
		prm_strucs_context.get_name_sets()
	);
}

/// \brief Get the domain corresponding to the specified index in the specified strucs_context
///
/// \pre `prm_index < size( prm_strucs_context )`
///
/// \relates strucs_context
domain_opt cath::file::get_domain_opt_of_index(const strucs_context &prm_strucs_context, ///< The strucs_context to query
                                               const size_t         &prm_index           ///< The index of the entry to query
                                               ) {
	const auto &name_set = prm_strucs_context.get_name_sets()[ prm_index ];
	const auto &regions  = prm_strucs_context.get_regions  ()[ prm_index ];

	return make_optional_if_fn(
		static_cast<bool>( regions ),
		[&] {
			const str_opt &domain_name_opt = name_set.get_domain_name_from_regions();
			return
				domain_name_opt
					? domain{ *regions, *domain_name_opt }
					: domain{ *regions                   };
		}
	);
}

/// \brief Hash the details of the specified strucs_context into the specified seed value
///
/// \relates strucs_context
void cath::file::non_crypto_hash(size_t               &prm_init_hash_value, ///< The initial hash seed
                                 const strucs_context &prm_strucs_context   /// The strucs_context to include in the hash
                                 ) {
	for (const size_t &index : indices( size( prm_strucs_context ) ) ) {
		const pdb                &pdb      = prm_strucs_context.get_pdbs     ()[ index ];
		const name_set           &name_set = prm_strucs_context.get_name_sets()[ index ];
		const region_vec_opt     &regions  = prm_strucs_context.get_regions  ()[ index ];

		// Follow libstdc++'s lead of attempting to make the value used for none/nullopt an "unusual" value
		static constexpr auto NULL_HASH_VAL = static_cast<size_t>( -3333 );

		non_crypto_hash   ( prm_init_hash_value, name_set );
		hash_value_combine( prm_init_hash_value, regions ? hash<string>{}( to_string( *regions )  ) : NULL_HASH_VAL );
		hash_value_combine( prm_init_hash_value,           hash<size_t>{}( pdb.get_num_residues() )                 );
	}
}

/// \brief Hash the details of the specified strucs_context into the specified seed value
///
/// \relates strucs_context
size_t cath::file::non_crypto_hash_copy(size_t                prm_init_hash_value, ///< The initial hash seed
                                        const strucs_context &prm_strucs_context   ///< The strucs_context to include in the hash
                                        ) {
	non_crypto_hash( prm_init_hash_value, prm_strucs_context );
	return prm_init_hash_value;
}
