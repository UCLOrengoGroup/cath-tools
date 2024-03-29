/// \file
/// \brief The superposition_context class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPERPOSITION_CONTEXT_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPERPOSITION_CONTEXT_HPP

#include <optional>

#include <boost/operators.hpp>

#include "cath/alignment/alignment.hpp"
#include "cath/file/strucs_context.hpp"
#include "cath/superposition/superposition.hpp"

// clang-format off
namespace cath::align { class alignment_context; }
namespace cath::file { class name_set_list; }
namespace cath::opts { class data_dirs_spec; }
namespace cath::sup { class superposition_content_spec; }
// clang-format on

namespace cath::sup {

	/// \brief Store a superposition along with the context of the PDBs and ids of the actual structures being superposed
	///
	/// ATM, this is little more than a tuple<pdb_list, str_vec, superposition> with nice names and extra functionality
	/// of optionally storing an alignment
	class superposition_context final : private boost::equality_comparable<superposition_context> {
	private:
		/// \brief The superposition itself
		superposition  the_superposition;

		/// \brief TODOCUMENT
		file::strucs_context context;

		/// \brief An optional alignment corresponding to the superposition
		::std::optional<align::alignment> any_alignment;

	public:
		superposition_context(superposition,
		                      file::strucs_context);

		superposition_context(superposition,
		                      file::strucs_context,
		                      align::alignment);

		superposition_context(superposition,
		                      file::pdb_list,
		                      file::name_set_list,
		                      chop::region_vec_opt_vec);

		superposition_context(superposition,
		                      file::pdb_list,
		                      file::name_set_list,
		                      chop::region_vec_opt_vec,
		                      align::alignment);

		[[nodiscard]] const superposition &       get_superposition() const;
		[[nodiscard]] const file::strucs_context &get_strucs_context() const;
		[[nodiscard]] bool                        has_alignment() const;
		[[nodiscard]] const align::alignment &    get_alignment() const;

		superposition_context & set_pdbs(const file::pdb_list &);
	};

	const file::pdb_list & get_pdbs(const superposition_context &);
	const file::name_set_list & get_name_sets(const superposition_context &);
	const chop::region_vec_opt_vec & get_regions(const superposition_context &);

	file::pdb_list get_restricted_pdbs(const superposition_context &);

	file::pdb get_supn_content_pdb(const file::pdb &,
	                               const chop::region_vec_opt &,
	                               const superposition_content_spec &);

	file::pdb_list get_supn_content_pdbs(const superposition_context &,
	                                     const superposition_content_spec &);

	superposition_context set_pdbs_copy(superposition_context,
	                                    const file::pdb_list &);

//	bool operator==(const superposition_context &,
//	                const superposition_context &);

//	std::ostream & operator<<(std::ostream &,
//	                          const superposition_context &);

	size_t get_num_entries(const superposition_context &);

	void load_pdbs_from_names(superposition_context &,
	                          const opts::data_dirs_spec &);

	superposition_context load_pdbs_from_names_copy(superposition_context,
	                                                const opts::data_dirs_spec &);

	align::alignment_context make_restricted_alignment_context(const superposition_context &);

	superposition_context superposition_context_from_ptree(const boost::property_tree::ptree &);

	void save_to_ptree(boost::property_tree::ptree &,
	                   const superposition_context &);

} // namespace cath::sup

namespace cath::common {

	/// \brief Specialisation of cath::common::read_from_ptree for superposition_context
	template <>
	inline sup::superposition_context read_from_ptree<sup::superposition_context>(const boost::property_tree::ptree &prm_ptree ///< The ptree from which to read the superposition_context
																					) {
		return sup::superposition_context_from_ptree( prm_ptree );
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPERPOSITION_CONTEXT_HPP
