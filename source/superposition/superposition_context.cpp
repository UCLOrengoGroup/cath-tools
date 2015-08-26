/// \file
/// \brief The superposition_context class definitions

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

#include "superposition_context.h"

#include <boost/log/trivial.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>

#include "alignment/alignment_context.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "superposition/io/superposition_io.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::sup;
using namespace cath::sup::detail;
using namespace std;

using boost::log::trivial::severity_level;
using boost::property_tree::json_parser::write_json;
using boost::property_tree::ptree;
using boost::range::for_each;

/// \brief TODOCUMENT
superposition_context::superposition_context(const pdb_list      &arg_pdbs,         ///< TODOCUMENT
                                             const str_vec       &arg_names,        ///< TODOCUMENT
                                             const superposition &arg_superposition ///< TODOCUMENT
                                             ) : pdbs(arg_pdbs),
                                                 names(arg_names),
                                                 the_superposition(arg_superposition) {
}


/// \brief TODOCUMENT
superposition_context::superposition_context(const pdb_list      &arg_pdbs,          ///< TODOCUMENT
                                             const str_vec       &arg_names,         ///< TODOCUMENT
                                             const superposition &arg_superposition, ///< TODOCUMENT
                                             const alignment     &arg_alignment      ///< TODOCUMENT
                                             ) : pdbs(arg_pdbs),
                                                 names(arg_names),
                                                 the_superposition(arg_superposition),
                                                 any_alignment(1, arg_alignment) {

}

/// \brief TODOCUMENT
const pdb_list & superposition_context::get_pdbs_cref() const {
	return pdbs;
}

/// \brief TODOCUMENT
const str_vec & superposition_context::get_names_cref() const {
	return names;
}

/// \brief TODOCUMENT
const superposition & superposition_context::get_superposition_cref() const {
	return the_superposition;
}

/// \brief TODOCUMENT
bool superposition_context::has_alignment() const {
	return static_cast<bool>( any_alignment );
}

/// \brief TODOCUMENT
const alignment & superposition_context::get_alignment_cref() const {
	if ( ! has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to get alignment from superposition_context that doesn't contain one"));
	}
	return *any_alignment;
}

/// \brief TODOCUMENT
///
/// \relates superposition_context
///
/// \relatesalso alignment_context
alignment_context cath::sup::make_alignment_context(const superposition_context &arg_superposition_context ///< TODOCUMENT
                                                    ) {
	if ( ! arg_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("WARNING: Unable to extract alignment from superposition_context with no alignment"));
	}
	return alignment_context(
		arg_superposition_context.get_pdbs_cref(),
		arg_superposition_context.get_names_cref(),
		arg_superposition_context.get_alignment_cref()
	);
}

/// \brief TODOCUMENT
///
/// At present, this stores the names and the superposition but does nothing
/// with the alignment or the PDBs
///
/// \relates superposition_context
void cath::sup::save_to_ptree(ptree                       &arg_ptree,      ///< TODOCUMENT
                              const superposition_context &arg_sup_context ///< TODOCUMENT
                              ) {
	if ( arg_sup_context.has_alignment() ) {
		BOOST_LOG_TRIVIAL( warning ) << "Whilst converting a superposition_context to JSON, its alignment will be ignored because that is not currently supported";
	}

	const auto supn_ptree          = make_ptree_of( arg_sup_context.get_superposition_cref() );
	const auto trans_ptrees        = supn_ptree.get_child( superposition_io_consts::TRANSFORMATIONS_KEY );

	const auto entries_key = superposition_io_consts::ENTRIES_KEY;
	arg_ptree.put_child( entries_key, ptree{} );
	auto &entries_ptree = arg_ptree.get_child( entries_key );

	for_each(
		arg_sup_context.get_names_cref(),
		trans_ptrees,
		[&] (const string &name, const pair<string, ptree> &trans_ptree) {
			ptree entry_ptree;
			entry_ptree.put      ( superposition_io_consts::NAME_KEY,           name               );
			entry_ptree.put_child( superposition_io_consts::TRANSFORMATION_KEY, trans_ptree.second );
			entries_ptree.push_back( make_pair( "", entry_ptree ) );
		}
	);
}

/// \brief TODOCUMENT
///
/// At present, this stores the names and the superposition but does nothing
/// with the alignment or the PDBs
///
/// \relates superposition_context
ptree cath::sup::make_ptree_of(const superposition_context &arg_sup_context ///< TODOCUMENT
                               ) {
	ptree new_ptree;
	save_to_ptree( new_ptree, arg_sup_context );
	return new_ptree;
}

/// \brief Create a JSON string to represent the specified superposition
///
/// At present, this stores the names and the superposition but does nothing
/// with the alignment or the PDBs
///
/// \relates superposition_context
string cath::sup::to_json_string(const superposition_context &arg_sup_context, ///< The superposition_context to represent in the JSON string
                                 const bool                  &arg_pretty_print ///< Whether to use whitespace (including line breaks) in the JSON to make it more human-readable
                                 ) {
	ostringstream json_ss;
	ptree temp_ptree;
	save_to_ptree( temp_ptree, arg_sup_context );
	write_json( json_ss, temp_ptree, arg_pretty_print );
	return json_ss.str();
}
