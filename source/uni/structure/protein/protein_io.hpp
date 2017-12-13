/// \file
/// \brief The protein_io class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_IO_HPP
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_IO_HPP

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "chopping/chopping_type_aliases.hpp"
#include "common/type_aliases.hpp"
#include "file/name_set/name_set.hpp"

#include <iostream>

namespace cath { class protein; }
namespace cath { namespace file { class dssp_file; } }
namespace cath { namespace file { class pdb; } }
namespace cath { namespace file { class sec_file; } }
namespace cath { namespace file { class wolf_file; } }
namespace cath { namespace file { enum class dssp_skip_policy : char; } }

namespace cath {

	protein read_protein_from_wolf_and_sec_files(const boost::filesystem::path &,
	                                             const boost::filesystem::path &,
	                                             const std::string & = "",
	                                             const ostream_ref_opt & = std::ref( std::cerr ) );

	protein read_protein_from_dssp_pdb_and_sec_files(const boost::filesystem::path &,
	                                                 const boost::filesystem::path &,
	                                                 const boost::filesystem::path &,
	                                                 const file::dssp_skip_policy &,
	                                                 const std::string & = "",
	                                                 const ostream_ref_opt & = std::ref( std::cerr ) );

	protein read_protein_from_dssp_and_pdb_and_calc_sec(const boost::filesystem::path &,
	                                                    const boost::filesystem::path &,
	                                                    const file::dssp_skip_policy &,
	                                                    const std::string & = "",
	                                                    const ostream_ref_opt & = std::ref( std::cerr ) );

	protein read_protein_from_dssp_and_pdb(const boost::filesystem::path &,
	                                       const boost::filesystem::path &,
	                                       const file::dssp_skip_policy &,
	                                       const std::string & = "",
	                                       const ostream_ref_opt & = boost::none );

	protein read_protein_from_pdb(const boost::filesystem::path &,
	                              const std::string & = "",
	                              const chop::region_vec_opt & = boost::none,
	                              std::ostream & = std::cerr);



	protein make_protein_from_pdb_and_calc_dssp(const file::pdb &,
	                                            const std::string &,
	                                            const ostream_ref_opt & = boost::none );

	protein make_protein_from_pdb_and_calc_dssp_and_sec(const file::pdb &,
	                                                    const std::string &,
	                                                    const ostream_ref_opt & = boost::none );

	protein read_protein_from_pdb_and_calc_dssp(const boost::filesystem::path &,
	                                            const std::string & = "",
	                                            const chop::region_vec_opt & = boost::none,
	                                            const ostream_ref_opt & = boost::none );

	protein read_protein_from_pdb_and_calc_dssp_and_sec(const boost::filesystem::path &,
	                                                    const std::string & = "",
	                                                    const chop::region_vec_opt & = boost::none,
	                                                    const ostream_ref_opt & = boost::none );



	protein make_protein_from_wolf_and_sec(const file::wolf_file &,
	                                       const file::sec_file &,
	                                       const std::string & = "",
	                                       const ostream_ref_opt & = std::ref( std::cerr ) );

	protein make_protein_from_dssp_pdb_and_sec(const file::dssp_file &,
	                                           const file::pdb &,
	                                           const file::sec_file &,
	                                           const file::dssp_skip_policy &,
	                                           const std::string & = "",
	                                           const ostream_ref_opt & = std::ref( std::cerr ) );

	void add_name_and_paint_sec_file_onto_protein(protein &,
	                                              const file::sec_file &,
	                                              const file::name_set & = file::name_set{},
	                                              const ostream_ref_opt & = std::ref( std::cerr ) );

	protein add_name_and_paint_sec_file_onto_protein_copy(protein,
	                                                      const file::sec_file &,
	                                                      const file::name_set & = file::name_set{},
	                                                      const ostream_ref_opt & = std::ref( std::cerr ) );

	void paint_sec_file_onto_protein(protein &,
	                                 const file::sec_file &,
	                                 const ostream_ref_opt &);

	protein paint_sec_file_onto_protein_copy(protein,
	                                         const file::sec_file &,
	                                         const ostream_ref_opt &);

	void calc_and_paint_sec_file_onto_protein(protein &,
	                                          const ostream_ref_opt &);

	protein calc_and_paint_sec_file_onto_protein_copy(protein,
	                                                  const ostream_ref_opt &);

	void remove_domin_res(protein &,
	                      const boost::filesystem::path &,
	                      const ostream_ref_opt &);

	void remove_domin_res(protein &,
	                      const size_size_pair_vec &,
	                      const ostream_ref_opt &);
} // namespace cath

#endif
