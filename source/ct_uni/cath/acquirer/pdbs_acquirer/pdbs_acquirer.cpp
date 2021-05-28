/// \file
/// \brief The pdbs_acquirer class definitions

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

#include "pdbs_acquirer.hpp"

#include <optional>

#include "cath/acquirer/pdbs_acquirer/file_list_pdbs_acquirer.hpp"
#include "cath/acquirer/pdbs_acquirer/istream_pdbs_acquirer.hpp"
#include "cath/chopping/domain/domain.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/clone/check_uptr_clone_against_this.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/logger.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_list.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/file/strucs_context.hpp"
#include "cath/options/options_block/pdb_input_options_block.hpp"
#include "cath/options/options_block/pdb_input_spec.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;

using ::std::cbegin;
using ::std::cend;
using ::std::cerr;
using ::std::istream;
using ::std::make_pair;
using ::std::nullopt;
using ::std::pair;
using ::std::unique_ptr;
using ::std::vector;

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<pdbs_acquirer> pdbs_acquirer::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief TODOCUMENT
///
/// \TODO Consider taking an ostream_ref_opt argument rather than assuming cerr
///       (fix all errors, *then* provide default of ::std::nullopt)
pdb_list_name_set_list_pair pdbs_acquirer::get_pdbs_and_names(istream    &prm_istream,                ///< TODOCUMENT
                                                              const bool &prm_remove_partial_residues ///< TODOCUMENT
                                                              ) const {
	pair<pdb_list, name_set_list> pdbs_and_names = do_get_pdbs_and_names( prm_istream );
	// Create a vector of PDBs to be superposed
	pdb_list      &pdbs  = pdbs_and_names.first;
	name_set_list &names = pdbs_and_names.second;

	// Check the number of source files and then grab them
//	if (pdbs.size() < 2) {
//		logger::log_and_exit(
//			logger::return_code::TOO_FEW_PDBS_FOR_ALIGNMENT,
//			"ERROR: There aren't enough PDBs in the input to perform this alignment/superposition"
//		);
//	}

	// If the number of names doesn't match the number of PDBs then throw a wobbly
	if ( names.size() != pdbs.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"The number of names ("
			+ std::to_string( names.size() )
			+ ") doesn't match the number of PDBs ("
			+ std::to_string( pdbs.size() )
			+ ")"
		));
	}

	return prm_remove_partial_residues ? make_pair( pdb_list_of_backbone_complete_subset_pdbs( pdbs, ref( cerr ) ), names )
	                                   : pdbs_and_names;
}

/// \brief Construct suitable pdbs_acquirer objects implied by the specified pdb_input_spec
///
/// NOTE: Keep this code in sync with get_num_acquirers()
///
/// \relates pdb_input_spec
uptr_vec<pdbs_acquirer> cath::opts::get_pdbs_acquirers(const pdb_input_spec &prm_pdb_input_spec ///< The pdb_input_spec to query
                                                       ) {
	vector<unique_ptr<pdbs_acquirer>> pdb_acquirers;
	if ( prm_pdb_input_spec.get_read_from_stdin() ) {
		pdb_acquirers.push_back( ::std::make_unique<istream_pdbs_acquirer>() );
	}
	if ( ! prm_pdb_input_spec.get_input_files().empty() ) {
		pdb_acquirers.push_back( ::std::make_unique<file_list_pdbs_acquirer>( prm_pdb_input_spec.get_input_files() ) );
	}

	if ( pdb_acquirers.size() != get_num_acquirers( prm_pdb_input_spec ) ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"The number of PDB acquirers "           + ::std::to_string( pdb_acquirers.size() )
			+ " doesn't match the expected number (" + ::std::to_string( get_num_acquirers( prm_pdb_input_spec ) ) + ")"
		));
	}

	return pdb_acquirers;
}

/// \brief Construct suitable pdbs_acquirer objects implied by the specified pdb_input_options_block
///
/// \relates pdb_input_options_block
uptr_vec<pdbs_acquirer> cath::opts::get_pdbs_acquirers(const pdb_input_options_block &prm_pdb_input_options_block ///< The pdb_input_options_block to query
                                                       ) {
	return get_pdbs_acquirers( prm_pdb_input_options_block.get_pdb_input_spec() );
}

/// \brief Construct the single pdbs_acquirer implied by the specified pdb_input_spec
///        (or throw an invalid_argument_exception if fewer/more are implied)
///
/// \relates pdb_input_spec
unique_ptr<pdbs_acquirer> cath::opts::get_pdbs_acquirer(const pdb_input_spec &prm_pdb_input_spec ///< The pdb_input_spec to query
                                                        ) {
	const auto pdbs_acquirers = get_pdbs_acquirers( prm_pdb_input_spec );
	if ( pdbs_acquirers.size() != 1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to get pdbs_acquirer failed because the number of pdb_acquirers isn't one"));
	}
	return pdbs_acquirers.front()->clone();
}

/// \brief Strip the specified domain_vec into a str_opt_vec and a region_vec_opt_vec,
///        making some attempt to avoid unnecessary allocations
///
/// There currently allocates the region_vec. More could be done to improve efficiency given motivation.
static pair<str_opt_vec, region_vec_opt_vec> strip_domain_vec(const domain_vec &prm_domains ///< The domain_vec to strip
                                                              ) {
	pair<str_opt_vec, region_vec_opt_vec> result;
	result.first .reserve( prm_domains.size() );
	result.second.reserve( prm_domains.size() );
	for (const domain &the_domain : prm_domains) {
		result.first .push_back( the_domain.get_opt_domain_id() );
		result.second.push_back( region_vec{ cbegin( the_domain ), cend( the_domain ) } );
	}
	return result;
}

/// \brief Combine the PDBs and names obtained from a pdbs_acquirer with IDs and regions to make a strucs_context
///
/// This is provided as a common place for this bit of behaviour to be performed so that it will
/// be easier to change in one place in the future
///
/// At present, this uses the prm_ids if they have the correct number of entries and prm_names_from_acq otherwise.
///
/// In the future, it may be worth building more interesting types (than str_vec) to record both the provenance (prm_names_from_acq)
/// and user-specified names (prm_ids) of the structure
strucs_context cath::opts::combine_acquired_pdbs_and_names_with_ids_and_domains(pdb_list          prm_pdbs,           ///< The PDBs obtained from a pdbs_acquirer
                                                                                name_set_list     prm_names_from_acq, ///< The names obtained from a pdbs_acquirer
                                                                                str_vec           specified_ids,      ///< Alternative IDs
                                                                                const domain_vec &prm_domains         ///< Regions for the strucs_context
                                                                                ) {
	auto stripped_domain_vec = strip_domain_vec( prm_domains );

	stripped_domain_vec.second.resize( prm_names_from_acq.size(), nullopt );

	return {
		std::move( prm_pdbs ),
		add_domain_names_from_regions_copy(
			add_specified_ids_copy(
				std::move( prm_names_from_acq ),
				std::move( specified_ids )
			),
			std::move( stripped_domain_vec.first )
		),
		std::move( stripped_domain_vec.second )
	};
}

/// \brief Use the specified PDBs acquirer to get PDBs and names and combine with the specified IDs and domains
///        to create a strucs_context
///
/// \relates pdbs_acquirer
strucs_context cath::opts::get_strucs_context(const pdbs_acquirer &prm_pdbs_acquirer,           ///< The pdbs_acquirer to use to get the PDBs
                                              istream             &prm_istream,                 ///< The istream, which may contain PDB data
                                              const bool          &prm_remove_partial_residues, ///< Whether to remove partial residues from the PDB data
                                              const str_vec       &prm_ids,                     ///< The IDs to set on the acquired strucs_context
                                              const domain_vec    &prm_domains                  ///< The domains to set on the acquired strucs_context
                                              ) {
	// Grab the PDBs and names
	auto pdbs_and_names = prm_pdbs_acquirer.get_pdbs_and_names(
		prm_istream,
		prm_remove_partial_residues
	);

	return combine_acquired_pdbs_and_names_with_ids_and_domains(
		std::move( pdbs_and_names.first  ),
		std::move( pdbs_and_names.second ),
		prm_ids,
		prm_domains
	);
}
