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

#include "display_colour_spec.h"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/join.hpp>

#include "alignment/alignment_context.h"
#include "common/algorithm/contains.h"
#include "common/algorithm/copy_build.h"
#include "common/algorithm/sort_uniq_build.h"
#include "common/size_t_literal.h"
#include "display/viewer/viewer.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_residue.h"
#include "superposition/superposition_context.h"

#include <algorithm>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::sup;
using namespace std;

using boost::adaptors::filtered;
using boost::adaptors::map_keys;
using boost::adaptors::map_values;
using boost::range::join;
using boost::lexical_cast;
using boost::none;
using boost::numeric_cast;

/// \brief TODOCUMENT
void display_colour_spec::colour_base(const display_colour &arg_colour,   ///< TODOCUMENT
                                      const bool           &arg_overwrite ///< TODOCUMENT
                                      ) {
	if (! arg_overwrite && base_clr ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to re-specify base colour"));
	}
	base_clr = arg_colour;
}

/// \brief TODOCUMENT
void display_colour_spec::colour_pdb(const size_t         &arg_pdb_index, ///< TODOCUMENT
		                             const display_colour &arg_colour,    ///< TODOCUMENT
		                             const bool           &arg_overwrite  ///< TODOCUMENT
		                             ) {
	if ( ! arg_overwrite && contains ( clr_of_pdb, arg_pdb_index ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to re-specify colour of PDB"));
	}
//	clr_of_pdb.insert( make_pair( arg_pdb_index, arg_colour ) );
	clr_of_pdb[ arg_pdb_index ] = arg_colour;
}

/// \brief TODOCUMENT
void display_colour_spec::colour_pdb_residue(const size_t         &arg_pdb_index,     ///< TODOCUMENT
                                             const size_t         &arg_residue_index, ///< TODOCUMENT
                                             const display_colour &arg_colour,        ///< TODOCUMENT
                                             const bool           &arg_overwrite      ///< TODOCUMENT
                                             ) {
	const size_size_pair pdb_and_residue( arg_pdb_index, arg_residue_index );
	if ( ! arg_overwrite && contains( clr_of_pdb_and_res, pdb_and_residue ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to re-specify colour of residue"));
	}
//	clr_of_pdb_and_res.insert( make_pair( pdb_and_residue, arg_colour ) );
	clr_of_pdb_and_res[ pdb_and_residue ] = arg_colour;
}

/// \brief TODOCUMENT
const opt_display_colour & display_colour_spec::get_base_clr() const {
	return base_clr;
}

/// \brief TODOCUMENT
const size_display_colour_map & display_colour_spec::get_clr_of_pdb() const {
	return clr_of_pdb;
}

/// \brief TODOCUMENT
const size_size_display_colour_map & display_colour_spec::get_clr_of_pdb_and_res() const {
	return clr_of_pdb_and_res;
}

/// \brief TODOCUMENT
///
/// relates display_colour_spec
opt_display_colour cath::get_clr_of_pdb_index(const display_colour_spec &arg_display_colour_spec, ///< TODOCUMENT
                                              const size_t              &arg_entry                ///< TODOCUMENT
                                              ) {
	const size_display_colour_map &colour_of_pdb = arg_display_colour_spec.get_clr_of_pdb();
	const bool                     has_entry     = contains( colour_of_pdb, arg_entry );
	return has_entry ? opt_display_colour( colour_of_pdb.at( arg_entry ) )
	                 : opt_display_colour( none                          );
}

/// \brief TODOCUMENT
///
/// relates display_colour_spec
opt_display_colour cath::get_clr_of_pdb_and_res_indices(const display_colour_spec &arg_display_colour_spec, ///< TODOCUMENT
                                                        const size_t              &arg_entry,               ///< TODOCUMENT
                                                        const size_t              &arg_index                ///< TODOCUMENT
                                                        ) {
	const size_size_pair                entry_index_pair = make_pair( arg_entry, arg_index );
	const size_size_display_colour_map &colour_of_pdb    = arg_display_colour_spec.get_clr_of_pdb_and_res();
	const bool                          has_entry        = contains( colour_of_pdb, entry_index_pair );
	return has_entry ? opt_display_colour( colour_of_pdb.at( entry_index_pair ) )
	                 : opt_display_colour( none                                 );
}

/// \brief TODOCUMENT
///
/// relates display_colour_spec
display_colour_vec cath::get_pdb_colours(const display_colour_spec &arg_colour_spec ///< TODOCUMENT
                                         ) {
	return sort_uniq_build<display_colour_vec>(
		arg_colour_spec.get_clr_of_pdb() | map_values
	);
}

/// \brief TODOCUMENT
///
/// relates display_colour_spec
display_colour_vec cath::get_residue_colours(const display_colour_spec &arg_colour_spec ///< TODOCUMENT
                                             ) {
	return sort_uniq_build<display_colour_vec>(
		arg_colour_spec.get_clr_of_pdb_and_res() | map_values
	);
}

/// \brief TODOCUMENT
///
/// relates display_colour_spec
display_colour_vec cath::get_all_colours(const display_colour_spec &arg_colour_spec ///< TODOCUMENT
                                         ) {
	const display_colour_vec pdb_colours     = get_pdb_colours    ( arg_colour_spec );
	const display_colour_vec residue_colours = get_residue_colours( arg_colour_spec );
	return sort_uniq_build<display_colour_vec>( join( pdb_colours, residue_colours ) );
}

/// \brief TODOCUMENT
///
/// relates display_colour_spec
size_vec cath::get_pdbs_of_colour(const display_colour_spec &arg_colour_spec, ///< TODOCUMENT
                                  const display_colour      &arg_colour       ///< TODOCUMENT
                                  ) {
	return copy_build<size_vec>(
		arg_colour_spec.get_clr_of_pdb()
		| filtered( [&] (const size_display_colour_map_val &x) { return ( x.second == arg_colour ); } )
		| map_keys
	);

}

/// \brief TODOCUMENT
///
/// relates display_colour_spec
size_size_vec_map cath::get_residues_of_colour(const display_colour_spec &arg_colour_spec, ///< TODOCUMENT
                                               const display_colour      &arg_colour       ///< TODOCUMENT
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
/// \brief TODOCUMENT
str_vec cath::generate_colour_names(const size_t &arg_num_colours
                                    ) {
	const size_t num_width = lexical_cast<string>( max( static_cast<size_t>( 1_z ), arg_num_colours ) - 1 ).length();
	str_vec colour_names;
	colour_names.reserve( arg_num_colours );
	for (size_t colour_ctr = 0; colour_ctr < arg_num_colours; ++colour_ctr) {
		ostringstream colour_name_ss;
		colour_name_ss << "cath_tools_defined_colour_";
		colour_name_ss << setw( numeric_cast<int>( num_width ) ) << setfill( '0' );
		colour_name_ss << colour_ctr;
		colour_names.push_back( colour_name_ss.str() );
	}
	return colour_names;
}

/// \brief TODOCUMENT
///
/// relates display_colour_spec
void cath::colour_viewer_with_spec(const display_colour_spec &arg_colour_spec,       ///< TODOCUMENT
                                   const viewer              &arg_viewer,            ///< TODOCUMENT
                                   const alignment_context   &arg_alignment_context, ///< TODOCUMENT
                                   ostream                   &arg_os                 ///< TODOCUMENT
                                   ) {
	const pdb_list          &pdbs          = arg_alignment_context.get_pdbs();
	const str_vec            cleaned_names = clean_names_for_viewer( arg_alignment_context );

	const display_colour_vec colours      = get_all_colours( arg_colour_spec );
	const size_t             num_colours  = colours.size();
	const str_vec            colour_names = generate_colour_names( num_colours );

	for (size_t colour_ctr = 0; colour_ctr < num_colours; ++colour_ctr) {
		arg_viewer.define_colour(
			arg_os,
			colours      [ colour_ctr ],
			colour_names [ colour_ctr ]
		);
	}

	for (size_t colour_ctr = 0; colour_ctr < num_colours; ++colour_ctr) {
		const display_colour &colour = colours[ colour_ctr ];
		const size_vec pdb_indices = get_pdbs_of_colour( arg_colour_spec, colour );
		for (const size_t &pdb_index : pdb_indices) {
			arg_viewer.colour_pdb(
				arg_os,
				colour_names [ colour_ctr ],
				cleaned_names[ pdb_index  ]
			);
		}
	}

	for (size_t colour_ctr = 0; colour_ctr < num_colours; ++colour_ctr) {
		const size_size_vec_map pdb_and_res_indices = get_residues_of_colour( arg_colour_spec, colours[ colour_ctr ] );
		for (const size_size_vec_pair &pdb_and_res_index : pdb_and_res_indices) {
			const size_t      &pdb_index       = pdb_and_res_index.first;
			const size_vec    &residue_indices = pdb_and_res_index.second;
			const pdb         &the_pdb         = pdbs[ pdb_index  ];

			residue_name_vec residue_names;
			residue_names.reserve( residue_indices.size() );
			for (const size_t &residue_index : residue_indices) {
				residue_names.push_back( the_pdb.get_residue_cref_of_backbone_complete_index( residue_index ).get_residue_name() ) ;
			}
			arg_viewer.colour_pdb_residues(
				arg_os,
				colour_names [ colour_ctr ],
				cleaned_names[ pdb_index  ],
				residue_names
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
