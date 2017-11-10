/// \file
/// \brief The display_colour_spec class definitions

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

#include "display_colour_spec.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/map.hpp>

#include <boost/range/join.hpp>

#include "alignment/alignment_context.hpp"
#include "chopping/region/region.hpp"
#include "common/algorithm/contains.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/sort_uniq_build.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "display/viewer/viewer.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_residue.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;

using boost::adaptors::filtered;
using boost::adaptors::map_keys;
using boost::adaptors::map_values;
using boost::none;
using boost::range::join;
using std::make_pair;
using std::ostream;
using std::ostringstream;
using std::string;

/// \brief Ctor from a broad_display_colour_spec
display_colour_spec::display_colour_spec(broad_display_colour_spec arg_broad_spec ///< The broad_display_colour_spec from which to construct this display_colour_spec
                                         ) : the_broad_spec { std::move( arg_broad_spec ) } {
}

/// \brief Specify the base colour to use
void display_colour_spec::colour_base(const display_colour &arg_colour,   ///< The colour to use
                                      const bool           &arg_overwrite ///< Whether to allow any previous base value to be overwritten without throwing
                                      ) {
	the_broad_spec.colour_base( arg_colour, arg_overwrite );
}

/// \brief Specify the colour to use for a specific structure
void display_colour_spec::colour_pdb(const size_t         &arg_struc_idx, ///< The index of the structure for which the colour should be set
		                             const display_colour &arg_colour,    ///< The colour to use
		                             const bool           &arg_overwrite  ///< Whether to allow any previous base value to be overwritten without throwing
		                             ) {
	the_broad_spec.colour_pdb( arg_struc_idx, arg_colour, arg_overwrite );
}

/// \brief Specify the colour to use for a specific residue
void display_colour_spec::colour_pdb_residue(const size_t         &arg_struc_idx,     ///< The index of the structure for which the colour should be set
                                             const size_t         &arg_residue_index, ///< The index of the residue within the structure for which the colour should be set
                                             const display_colour &arg_colour,        ///< The colour to use
                                             const bool           &arg_overwrite      ///< Whether to allow any previous base value to be overwritten without throwing
                                             ) {
	const size_size_pair pdb_and_residue( arg_struc_idx, arg_residue_index );
	if ( ! arg_overwrite && contains( clr_of_pdb_and_res, pdb_and_residue ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to re-specify colour of residue"));
	}
	clr_of_pdb_and_res[ pdb_and_residue ] = arg_colour;
}

/// \brief Getter for the broad_display_colour_spec
const broad_display_colour_spec & display_colour_spec::get_broad_spec() const {
	return the_broad_spec;
}

/// \brief Getter for the colours of structure/residue
const size_size_display_colour_map & display_colour_spec::get_clr_of_pdb_and_res() const {
	return clr_of_pdb_and_res;
}

/// \brief Getter for the base colour
///
/// relates display_colour_spec
const display_colour_opt & cath::get_base_clr(const display_colour_spec &arg_display_colour_spec ///< The display_colour_spec to query
                                              ) {
	return arg_display_colour_spec.get_broad_spec().get_base_clr();
}

/// \brief Getter for colour by structure index
///
/// relates display_colour_spec
const size_display_colour_map & cath::get_clr_of_pdb(const display_colour_spec &arg_display_colour_spec ///< The display_colour_spec to query
                                                     ) {
	return arg_display_colour_spec.get_broad_spec().get_clr_of_pdb();
}

/// \brief Get the (optional) colour for the structure of the specified index
///
/// relates display_colour_spec
display_colour_opt cath::get_clr_of_pdb_index(const display_colour_spec &arg_display_colour_spec, ///< The display_colour_spec to query
                                              const size_t              &arg_entry                ///< The index of the structure required
                                              ) {
	return get_clr_of_pdb_index( arg_display_colour_spec.get_broad_spec(), arg_entry );
}

/// \brief Get the colour of structure/residue of the specified indices
///
/// relates display_colour_spec
display_colour_opt cath::get_clr_of_pdb_and_res_indices(const display_colour_spec &arg_display_colour_spec, ///< The display_colour_spec to query
                                                        const size_t              &arg_entry,               ///< The index of the structure required
                                                        const size_t              &arg_index                ///< The index of the residue required within the structure
                                                        ) {
	const size_size_pair                entry_index_pair = make_pair( arg_entry, arg_index );
	const size_size_display_colour_map &colour_of_pdb    = arg_display_colour_spec.get_clr_of_pdb_and_res();
	const bool                          has_entry        = contains( colour_of_pdb, entry_index_pair );
	return has_entry ? display_colour_opt( colour_of_pdb.at( entry_index_pair ) )
	                 : display_colour_opt( none                                 );
}

/// \brief Get a list of all colours used to colour whole structures
///
/// relates display_colour_spec
display_colour_vec cath::get_pdb_colours(const display_colour_spec &arg_colour_spec ///< The display_colour_spec to query
                                         ) {
	return get_pdb_colours( arg_colour_spec.get_broad_spec() );
}

/// \brief Get list of all colours used to colour residues
///
/// relates display_colour_spec
display_colour_vec cath::get_residue_colours(const display_colour_spec &arg_colour_spec ///< The display_colour_spec to query
                                             ) {
	return sort_uniq_build<display_colour_vec>(
		arg_colour_spec.get_clr_of_pdb_and_res() | map_values
	);
}

/// \brief Get a list of all colours
///
/// relates display_colour_spec
display_colour_vec cath::get_all_colours(const display_colour_spec &arg_colour_spec ///< The display_colour_spec to query
                                         ) {
	const display_colour_vec pdb_colours     = get_pdb_colours    ( arg_colour_spec );
	const display_colour_vec residue_colours = get_residue_colours( arg_colour_spec );
	return sort_uniq_build<display_colour_vec>( boost::range::join( pdb_colours, residue_colours ) );
}

/// \brief Get a list of the indices of all PDBs of the specified colour
///
/// relates display_colour_spec
size_vec cath::get_pdbs_of_colour(const display_colour_spec &arg_colour_spec, ///< The display_colour_spec to query
                                  const display_colour      &arg_colour       ///< The colour for which to search
                                  ) {
	return get_pdbs_of_colour( arg_colour_spec.get_broad_spec(), arg_colour );
}

/// \brief Get a list of all residues of the specified colour
///
/// relates display_colour_spec
size_size_vec_map cath::get_residues_of_colour(const display_colour_spec &arg_colour_spec, ///< The display_colour_spec to query
                                               const display_colour      &arg_colour       ///< The colour for which to search
                                               ) {
	// Get a vector of all pdb/residue indices for all residues relating to the specified colour
	const auto separate_residues = copy_build<size_size_pair_vec>(
		arg_colour_spec.get_clr_of_pdb_and_res()
		| filtered( [&] (const size_size_display_colour_map_val &x) { return ( x.second == arg_colour ); } )
		| map_keys
	);

	// Build and return a map containing a vector of residues indices for each pdb index
	size_size_vec_map residues;
	for (const size_size_pair &res : separate_residues) {
		residues[ res.first ].push_back( res.second );
	}
	return residues;
}

/// \brief Write (possibly residue-specific) colouring instructions for the specified viewer to the specified ostream
///        using the specified display_colour_spec and (viewer-cleaned) names
///
/// relates display_colour_spec
void cath::colour_viewer_with_spec(const display_colour_spec &arg_colour_spec,       ///< The specification of colours to use
                                   const viewer              &arg_viewer,            ///< The viewer to specify how to render the instructions
                                   const alignment_context   &arg_alignment_context, ///< The (viewer-cleaned) names for the structures to be coloured
                                   ostream                   &arg_os                 ///< The ostream to which the instructions should be written
                                   ) {
	const display_colour_vec colours       = get_all_colours( arg_colour_spec );
	const str_vec            cleaned_names = clean_names_for_viewer( arg_alignment_context );

	detail::colour_base_and_pdbs_impl(
		colours,
		arg_colour_spec.get_broad_spec(),
		arg_viewer,
		get_pdbs( arg_alignment_context ),
		cleaned_names,
		colour_category::STRUC_OR_RES,
		arg_os
	);

	const pdb_list &pdbs         = get_pdbs( arg_alignment_context );
	const size_t    num_colours  = colours.size();

	for (const size_t &colour_ctr : indices( num_colours ) ) {
		const auto pdb_and_res_indices = get_residues_of_colour( arg_colour_spec, colours[ colour_ctr ] );
		for (const auto &pdb_and_res_index : pdb_and_res_indices) {
			const size_t &pdb_index = pdb_and_res_index.first;
			arg_os << arg_viewer.get_colour_pdb_residues_str(
				generate_colour_name( colour_ctr, num_colours,colour_category::STRUC_OR_RES ),
				cleaned_names[ pdb_index ],
				transform_build<residue_id_vec>(
					pdb_and_res_index.second,
					[&] (const size_t &x) {
						return get_residue_of_backbone_complete_index( pdbs[ pdb_index ], x ).get_residue_id();
					}
				)
			);
		}
	}
}

///// \brief TODOCUMENT
/////
///// relates display_colour_spec
//void cath::colour_alignment_with_spec(const display_colour_spec &arg_colour_spec, ///<
//	                                  const alignment           &arg_alignment,
//	                                  const pdb_list            &arg_pbs,
//	                                  const str_vec             &arg_names,
//	                                  ostream                   &arg_ostream
//	                                  ) {
//	"<span class=value3>S</span>";
//	"<span class=value3>S</span>";
//}
