/// \file
/// \brief The protein_source_file_set class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET_PROTEIN_SOURCE_FILE_SET_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET_PROTEIN_SOURCE_FILE_SET_HPP

#include <filesystem>
#include <optional>

#include <boost/ptr_container/ptr_vector.hpp>

#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/file/file_type_aliases.hpp"
#include "cath/structure/protein/protein_list.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_file_combn.hpp"

// clang-format off
namespace cath::opts { class data_dirs_spec; }
// clang-format on

namespace cath {

	/// \brief Represent a group of files from which to populate a protein (and the method of doing it)
	///
	/// The concrete classes need to define which files they require (do_get_file_set())
	/// as well as how to read those files into a protein (do_read_files()).
	///
	/// This is used by protein_source_file_set::read_files() to populate a data_file_path_map with the required filenames
	/// before passing that to do_read_files().
	class protein_source_file_set {
	private:
		/// \brief Pure virtual method with which each concrete protein_source_file_set must define how to create a clone of itself
		[[nodiscard]] virtual std::unique_ptr<protein_source_file_set> do_clone() const = 0;

		/// \brief Pure virtual method with which each concrete protein_source_file_set must define the list of files it reads
		[[nodiscard]] virtual file::data_file_vec do_get_file_set() const = 0;

		/// \brief Pure virtual method with which each concrete protein_source_file_set must define the primary file it reads (which must be one of those returned by do_get_file_set())
		[[nodiscard]] virtual file::data_file do_get_primary_file() const = 0;

		/// \brief Pure virtual method with which each concrete protein_source_file_set must define the equivalent protein_file_combn value
		[[nodiscard]] virtual protein_file_combn do_get_protein_file_combn() const = 0;

		/// \brief Pure virtual method with which each concrete protein_source_file_set must define whether it
		///        makes proteins that are SSAP-ready (with data loaded for sec, phi/psi accessibility etc)
		[[nodiscard]] virtual bool do_makes_ssap_ready_protein() const = 0;

		/// \brief Pure virtual method with which each concrete protein_source_file_set must define the behaviour to read the specified files
		///        from a data_file_path_map that is pre-populated with the file types the concrete protein_source_file_set requires
		virtual protein do_read_files(const file::data_file_path_map &,
		                              const std::string &,
		                              std::ostream &) const = 0;

		/// \brief Virtual method with which each concrete protein_source_file_set may define how to do_read_files() and
		///        restrict the result to the specified regions.
		///
		/// The default implementation uses restrict_to_regions_copy() but a concrete protein_source_file_set
		/// may override (and then probably implement do_read_files() terms of that)
		virtual protein do_read_and_restrict_files(const file::data_file_path_map &,
		                                           const std::string &,
		                                           const chop::region_vec_opt &,
		                                           std::ostream &) const;

	public:
		protein_source_file_set() = default;
		virtual ~protein_source_file_set() noexcept = default;

		protein_source_file_set(const protein_source_file_set &) = default;
		protein_source_file_set(protein_source_file_set &&) noexcept = default;
		protein_source_file_set & operator=(const protein_source_file_set &) = default;
		protein_source_file_set & operator=(protein_source_file_set &&) noexcept = default;
		
		[[nodiscard]] std::unique_ptr<protein_source_file_set> clone() const;

		[[nodiscard]] file::data_file_vec get_file_set() const;
		[[nodiscard]] file::data_file     get_primary_file() const;
		[[nodiscard]] protein_file_combn  get_protein_file_combn() const;
		[[nodiscard]] bool                makes_ssap_ready_protein() const;
		protein read_files( const opts::data_dirs_spec &, const std::string &, const chop::region_vec_opt &, std::ostream & ) const;
	};

	protein read_protein_from_files(const protein_source_file_set &,
	                                const opts::data_dirs_spec &,
	                                const std::string &,
	                                const chop::domain_opt &,
	                                std::ostream &);

	protein read_protein_from_files(const protein_source_file_set &,
	                                const ::std::filesystem::path &,
	                                const std::string &,
	                                const ostream_ref_opt & = ::std::nullopt);

	protein_list read_proteins_from_files(const protein_source_file_set &,
	                                      const ::std::filesystem::path &,
	                                      const str_vec &,
	                                      const ostream_ref_opt & = ::std::nullopt);

	file::data_file_path_map get_filename_of_data_file(const protein_source_file_set &,
	                                                   const opts::data_dirs_spec &,
	                                                   const std::string &);

	::std::filesystem::path get_primary_file_from_map(const protein_source_file_set &,
	                                                  const file::data_file_path_map &);

	using protein_source_file_set_pvec = boost::ptr_vector<protein_source_file_set>;
	using protein_file_combn_vec = std::vector<protein_file_combn>;

	protein_source_file_set_pvec get_all_protein_source_file_sets();
	protein_file_combn_vec get_ssap_ready_protein_file_combns(const protein_source_file_set_pvec &);
	protein_file_combn_vec get_protein_file_combns(const protein_source_file_set_pvec &);
	protein_file_combn_vec get_all_protein_file_combns();
	protein_file_combn_vec get_ssap_ready_protein_file_combns();
	str_vec get_protein_file_combn_strings(const protein_file_combn_vec &);
	str_vec get_ssap_ready_protein_file_combn_strings();
	str_vec get_all_protein_file_combn_strings();

 	/// \brief Function to make protein_source_file_set meet the Clonable concept (used in ptr_container)
 	///
 	/// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
 	///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
 	///
 	/// This gets the smart pointer from the clone() method and then calls release on it.
 	///
 	/// \returns A raw pointer to a new copy of the protein_source_file_set argument, with the same dynamic type.
 	///          The caller is responsible for deleting this new object.
 	inline protein_source_file_set * new_clone(const protein_source_file_set &prm_protein_source_file_set ///< The protein_source_file_set to clone
 	                                           ) {
 	        return prm_protein_source_file_set.clone().release();
 	}
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET_PROTEIN_SOURCE_FILE_SET_HPP
