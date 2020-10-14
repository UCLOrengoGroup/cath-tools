/// \file
/// \brief The alignment_context class definitions

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

#include "alignment_context.hpp"

#include "cath/acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "cath/acquirer/pdbs_acquirer/pdbs_acquirer.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/superposition/superposition_context.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::chop;
using namespace ::cath::file;
using namespace ::cath::opts;
using namespace ::cath::sup;

using ::std::istream;

/// \brief Ctor for alignment_context
alignment_context::alignment_context(alignment      prm_alignment, ///< TODOCUMENT
                                     strucs_context prm_context    ///< TODOCUMENT
                                     ) : the_alignment { std::move( prm_alignment ) },
                                         context       { std::move( prm_context   ) } {
}

/// \brief Ctor for alignment_context
alignment_context::alignment_context(alignment                 prm_alignment, ///< TODOCUMENT
                                     const pdb_list           &prm_pdbs,      ///< TODOCUMENT
                                     const name_set_list      &prm_name_sets, ///< TODOCUMENT
                                     const region_vec_opt_vec &prm_regions    ///< The specification of the regions of the PDBs to which the alignment refers
                                     ) : the_alignment { std::move( prm_alignment )       },
                                         context       { prm_pdbs, prm_name_sets, prm_regions } {
}

/// \brief TODOCUMENT
const alignment & alignment_context::get_alignment() const {
	return the_alignment;
}

/// \brief TODOCUMENT
const strucs_context & alignment_context::get_strucs_context() const {
	return context;
}

/// \brief TODOCUMENT
const pdb_list & cath::align::get_pdbs(const alignment_context &prm_align_context ///< TODOCUMENT
                                       ) {
	return prm_align_context.get_strucs_context().get_pdbs();
}

/// \brief TODOCUMENT
const name_set_list & cath::align::get_name_sets(const alignment_context &prm_align_context ///< TODOCUMENT
                                                 ) {
	return prm_align_context.get_strucs_context().get_name_sets();
}

/// \brief Getter for the specification of the regions of the PDBs to which the alignment refers
const region_vec_opt_vec & cath::align::get_regions(const alignment_context &prm_align_context ///< TODOCUMENT
                                                    ) {
	return prm_align_context.get_strucs_context().get_regions();
}

/// \brief Get the number of entries represented in the specified alignment_context
///
/// \relates alignment_context
size_t cath::align::get_num_entries(const alignment_context &prm_alignment_context ///< The alignment_context to query
                                    ) {
	return size( prm_alignment_context.get_strucs_context() );
}

/// \brief Get a copy of the PDBs in the specified alignment_context, restricted to its regions
///
/// \relates alignment_context
pdb_list cath::align::get_restricted_pdbs(const alignment_context &prm_alignment_context ///< The alignment_context to query
                                          ) {
  return get_restricted_pdbs( prm_alignment_context.get_strucs_context() );
}

/// \brief Make an alignment context of the specified alignment and strucs_context, with the
///        PDBs restricted according to the regions in the strucs_context
///
/// \relates alignment_context
alignment_context cath::align::make_restricted_alignment_context(alignment      prm_alignment,     ///< The alignment from which to build the alignment_context
                                                                 strucs_context prm_strucs_context ///< The strucs_context from which to build the alignment_context
                                                                 ) {
	return {
		std::move( prm_alignment ),
		restrict_pdbs_copy( std::move( prm_strucs_context ) )
	};
}

/// \brief TODOCUMENT
///
/// \relates alignment_context
///
/// \relatesalso superposition_context
superposition_context cath::align::make_superposition_context(const alignment_context &prm_alignment_context, ///< TODOCUMENT
                                                              const superposition     &prm_superposition      ///< TODOCUMENT
                                                              ) {
	return superposition_context(
		prm_superposition,
		prm_alignment_context.get_strucs_context(),
		prm_alignment_context.get_alignment()
	);
}

/// \brief Get an alignment context
///
/// \relates alignment_acquirer
alignment_context cath::align::get_alignment_context(const alignment_acquirer &prm_alignment_acquirer,      ///< The alignment_acquirer to use to get the alignment
                                                     const pdbs_acquirer      &prm_pdbs_acquirer,           ///< The pdbs_acquirer to use to get the PDBs
                                                     istream                  &prm_istream,                 ///< The istream, which may contain PDB data
                                                     const bool               &prm_remove_partial_residues, ///< Whether to remove partial residues from the PDB data
                                                     const str_vec            &prm_ids,                     ///< The IDs to set on the acquired strucs_context
                                                     const domain_vec         &prm_domains,                 ///< The domains to set on the acquired strucs_context
                                                     const align_refining     &prm_align_refining           ///< How much refining should be done to the alignment
                                                     ) {
	const auto the_strucs_context = get_strucs_context(
		prm_pdbs_acquirer,
		prm_istream,
		prm_remove_partial_residues,
		prm_ids,
		prm_domains
	);
	return {
		prm_alignment_acquirer
			.get_alignment_and_spanning_tree( restrict_pdbs_copy( the_strucs_context ), prm_align_refining )
			.first,
		the_strucs_context
	};
}
