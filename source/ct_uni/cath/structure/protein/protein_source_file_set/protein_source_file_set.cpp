/// \file
/// \brief The protein_source_file_set class definitions

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

#include "protein_source_file_set.hpp"

#include <filesystem>
#include <map>

#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/transform.hpp>

#include "cath/chopping/domain/domain.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/clone/check_uptr_clone_against_this.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/file/options/data_dirs_options_block.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_pdb.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_pdb_and_calc.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_pdb_and_dssp_and_calc.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_pdb_dssp_and_sec.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_wolf_and_sec.hpp"
#include "cath/structure/protein/residue.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;

using ::boost::adaptors::filtered;
using ::boost::assign::ptr_push_back;
using ::boost::lexical_cast;
using ::boost::none;
using ::boost::range::transform;
using ::std::back_inserter;
using ::std::filesystem::path;
using ::std::make_pair;
using ::std::ostream;
using ::std::ostringstream;
using ::std::string;
using ::std::unique_ptr;

/// \brief Read a PDB and restrict it by the specified regions
protein protein_source_file_set::do_read_and_restrict_files(const data_file_path_map &prm_filename_of_data_file, ///< The pre-loaded map of file types to filenames
                                                            const string             &prm_protein_name,          ///< The name of the protein that is to be read from files
                                                            const region_vec_opt     &prm_regions,               ///< The regions to which the resulting protein should be restricted
                                                            ostream                  &prm_stderr                 ///< The ostream to which any warnings/errors should be written
                                                            ) const {
	protein the_protein = do_read_files( prm_filename_of_data_file, prm_protein_name, prm_stderr );
	return prm_regions
		? restrict_to_regions_copy(
			std::move( the_protein ),
			prm_regions
		)
		: the_protein;
}

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<protein_source_file_set> protein_source_file_set::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief An NVI pass-through method to get the list of files to be read
data_file_vec protein_source_file_set::get_file_set() const {
	return do_get_file_set();
}

/// \brief An NVI pass-through method to get the primary file
data_file protein_source_file_set::get_primary_file() const {
	return do_get_primary_file();
}

/// \brief An NVI pass-through method to get the equivalent protein_file_combn value
protein_file_combn protein_source_file_set::get_protein_file_combn() const {
	return do_get_protein_file_combn();
}

/// \brief An NVI pass-through method to return whether this makes proteins that are SSAP-ready
///        (with data loaded for sec, phi/psi accessibility etc)
bool protein_source_file_set::makes_ssap_ready_protein() const {
	return do_makes_ssap_ready_protein();
}

/// \brief An NVI pass-through method to read the files that have been specified
protein protein_source_file_set::read_files(const data_dirs_spec &prm_data_dirs,    ///< The data_dirs_options_block to specify how things should be done
                                            const string         &prm_protein_name, ///< The name of the protein that is to be read from files
                                            const region_vec_opt &prm_regions,      ///< The regions to which the resulting protein should be restricted
                                            ostream              &prm_stderr        ///< The ostream to which any warnings/errors should be written
                                            ) const {
	const data_file_path_map filename_of_data_file = get_filename_of_data_file(
		*this,
		prm_data_dirs,
		prm_protein_name
	);

	return do_read_and_restrict_files(
		filename_of_data_file,
		prm_protein_name,
		prm_regions,
		prm_stderr
	).set_name_set( name_set{
		get_primary_file_from_map( *this, filename_of_data_file ),
		prm_protein_name
	} );
}

/// \brief Read a protein from the specified files and apply any name from the optional domain
///
/// \relates protein_source_file_set
protein cath::read_protein_from_files(const protein_source_file_set &prm_source_file_set, /// The protein_source_file_set specifying which set of files should be used to build the protein
                                      const data_dirs_spec          &prm_data_dirs,       ///< The data_dirs_options_block to specify how things should be done
                                      const string                  &prm_protein_name,    ///< The name of the protein that is to be read from files
                                      const domain_opt              &prm_domain,          ///< The domain to which the resulting protein should be restricted
                                      ostream                       &prm_stderr           ///< The ostream to which any warnings/errors should be written
                                      ) {
	const region_vec_opt regions = make_optional_if_fn(
		static_cast<bool>( prm_domain ),
		[&] { return region_vec{ common::cbegin( *prm_domain ), common::cend( *prm_domain ) }; }
	);
	protein the_protein = prm_source_file_set.read_files(
		prm_data_dirs,
		prm_protein_name,
		regions,
		prm_stderr
	);
	if ( prm_domain && has_domain_id( *prm_domain ) ) {
		the_protein.get_name_set().set_domain_name_from_regions( get_domain_id( *prm_domain ) );
	}
	return the_protein;
}

/// \brief TODOCUMENT
///
/// \relates protein_source_file_set
protein cath::read_protein_from_files(const protein_source_file_set &prm_source_file_set, /// The protein_source_file_set specifying which set of files should be used to build the protein
                                      const path                    &prm_data_dir,        ///< The directory from which the files should be read
                                      const string                  &prm_protein_name,    ///< The name of the protein that is to be read from files
                                      const ostream_ref_opt         &prm_ostream          ///< An optional reference to an ostream to which any warnings/errors should be written
                                      ) {
	ostringstream parse_ss;
	return prm_source_file_set.read_files(
		build_data_dirs_spec_of_dir( prm_data_dir ),
		prm_protein_name,
		none,
		( prm_ostream ? prm_ostream->get() : parse_ss )
	);
}

/// \brief TODOCUMENT
///
/// \relates protein_source_file_set
protein_list cath::read_proteins_from_files(const protein_source_file_set &prm_source_file_set, /// The protein_source_file_set specifying which set of files should be used to build the protein
                                            const path                    &prm_data_dir,        ///< The directory from which the files should be read
                                            const str_vec                 &prm_protein_names,   ///< The name of the protein that is to be read from files
                                            const ostream_ref_opt         &prm_ostream          ///< An optional reference to an ostream to which any warnings/errors should be written
                                            ) {
	protein_list the_proteins;
	transform(
		prm_protein_names,
		back_inserter( the_proteins ),
		[&] (const string &x) {
			return read_protein_from_files(
				prm_source_file_set,
				prm_data_dir,
				x,
				prm_ostream
			);
		}
	);
	return the_proteins;
}


/// \brief For each file type required by the specified protein_source_file_set, get the filenames associated with
///        the specified protein name in the specified data_dirs_spec
///
/// \returns A map of data_file to the filename selected for that file type
///
/// \relates protein_source_file_set
data_file_path_map cath::get_filename_of_data_file(const protein_source_file_set &prm_protein_source_file_set, ///< The protein_source_file_set which specifies the file types it requires for reading a protein
                                                   const data_dirs_spec          &prm_data_dirs,               ///< The data_dirs_spec specifying how to map a file type and a name to a filename
                                                   const string                  &prm_protein_name             ///< The name of the protein that is to be read from files
                                                   ) {
	const auto file_set = prm_protein_source_file_set.get_file_set();
	return transform_build<data_file_path_map>(
		file_set,
		[&] (const data_file &x) {
			return make_pair( x, find_file( prm_data_dirs, x, prm_protein_name ) );
		}
	);
}

/// \brief Return the file in the specified data_file_path_map corresponding to
///        the specified protein_source_file_set's primary file
///
/// \relates protein_source_file_set
path cath::get_primary_file_from_map(const protein_source_file_set &prm_protein_source_file_set, ///< The protein_source_file_set to query
                                     const data_file_path_map      &prm_filename_of_data_file    ///< The data_file_path_map to query
                                     ) {
	return prm_filename_of_data_file.at(
		prm_protein_source_file_set.get_primary_file()
	);
}

/// \brief Get all types of protein_source_file_set that are currently available
///
/// \relates protein_source_file_set
protein_source_file_set_pvec cath::get_all_protein_source_file_sets() {
	protein_source_file_set_pvec file_sets;
	ptr_push_back< protein_from_pdb                   >( file_sets )();
	ptr_push_back< protein_from_pdb_and_calc          >( file_sets )();
	ptr_push_back< protein_from_pdb_and_dssp_and_calc >( file_sets )();
	ptr_push_back< protein_from_pdb_dssp_and_sec      >( file_sets )();
	ptr_push_back< protein_from_wolf_and_sec          >( file_sets )();
	return file_sets;
}

/// \brief Convert a ptr_vector of protein_source_file_sets to a vector<protein_file_combn>
///
/// \relates protein_source_file_set
protein_file_combn_vec cath::get_ssap_ready_protein_file_combns(const protein_source_file_set_pvec &prm_protein_source_file_sets ///< The ptr_vector of protein_source_file_sets to convert
                                                                ) {
	return transform_build<protein_file_combn_vec>(
		prm_protein_source_file_sets
			| filtered( [] (const auto &x) { return x.makes_ssap_ready_protein(); } ),
		[] (const protein_source_file_set &x) { return x.get_protein_file_combn(); }
	);
}

/// \brief Convert a ptr_vector of protein_source_file_sets to a vector<protein_file_combn>
///
/// \relates protein_source_file_set
protein_file_combn_vec cath::get_protein_file_combns(const protein_source_file_set_pvec &prm_protein_source_file_sets ///< The ptr_vector of protein_source_file_sets to convert
                                                     ) {
	return transform_build<protein_file_combn_vec>(
		prm_protein_source_file_sets,
		[] (const protein_source_file_set &x) { return x.get_protein_file_combn(); }
	);
}

/// \brief Get all types of protein_source_file_set, converted to a vector<protein_file_combn>
///
/// \relates protein_source_file_set
protein_file_combn_vec cath::get_ssap_ready_protein_file_combns() {
	return get_ssap_ready_protein_file_combns( get_all_protein_source_file_sets() );
}


/// \brief Get all types of protein_source_file_set, converted to a vector<protein_file_combn>
///
/// \relates protein_source_file_set
protein_file_combn_vec cath::get_all_protein_file_combns() {
	return get_protein_file_combns( get_all_protein_source_file_sets() );
}

/// \brief Lexically-cast the members of a vector<protein_file_combn> to get a str_vec
///
/// \todo Create a "lexically_casted" Range adaptor, which can do this more generally
///
/// \relates protein_source_file_set
str_vec cath::get_protein_file_combn_strings(const protein_file_combn_vec &prm_protein_file_combns ///< The protein_file_combn values to convert
                                             ) {
	return transform_build<str_vec>(
		prm_protein_file_combns,
		[] (const protein_file_combn &x) { return lexical_cast<string>( x ); }
	);
}

/// \brief Get all types of protein_source_file_set, converted to a vector<protein_file_combn> and then to a str_vec
///
/// \relates protein_source_file_set
str_vec cath::get_ssap_ready_protein_file_combn_strings() {
	return get_protein_file_combn_strings( get_ssap_ready_protein_file_combns() );
}

/// \brief Get all types of protein_source_file_set, converted to a vector<protein_file_combn> and then to a str_vec
///
/// \relates protein_source_file_set
str_vec cath::get_all_protein_file_combn_strings() {
	return get_protein_file_combn_strings( get_all_protein_file_combns() );
}

