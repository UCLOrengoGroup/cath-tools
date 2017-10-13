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

#ifndef _CATH_TOOLS_SOURCE_SUPERPOSITION_SUPERPOSITION_CONTEXT_H
#define _CATH_TOOLS_SOURCE_SUPERPOSITION_SUPERPOSITION_CONTEXT_H

#include <boost/operators.hpp>
#include <boost/optional.hpp>

#include "alignment/alignment.hpp"
#include "file/strucs_context.hpp"
#include "superposition/superposition.hpp"


namespace cath { namespace align { class alignment_context; } }
namespace cath { namespace file { class name_set_list; } }
namespace cath { namespace opts { class data_dirs_spec; } }
namespace cath { namespace sup { class superposition_content_spec; } }

namespace cath {
	namespace sup {

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
			boost::optional<align::alignment> any_alignment;

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

			const superposition & get_superposition() const;
			const file::strucs_context & get_strucs_context() const;
			bool has_alignment() const;
			const align::alignment & get_alignment() const;

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

//		bool operator==(const superposition_context &,
//		                const superposition_context &);

//		std::ostream & operator<<(std::ostream &,
//		                          const superposition_context &);

		size_t get_num_entries(const superposition_context &);

		void load_pdbs_from_names(superposition_context &,
		                          const opts::data_dirs_spec &);

		superposition_context load_pdbs_from_names_copy(superposition_context,
		                                                const opts::data_dirs_spec &);

		align::alignment_context make_restricted_alignment_context(const superposition_context &);

		superposition_context superposition_context_from_ptree(const boost::property_tree::ptree &);

		void save_to_ptree(boost::property_tree::ptree &,
		                   const superposition_context &);

	} // namespace sup

	namespace common {

		/// \brief Specialisation of cath::common::read_from_ptree for superposition_context
		template <>
		inline sup::superposition_context read_from_ptree<sup::superposition_context>(const boost::property_tree::ptree &arg_ptree ///< The ptree from which to read the superposition_context
		                                                                              ) {
			return sup::superposition_context_from_ptree( arg_ptree );
		}

	} // namespace common

} // namespace cath

#endif
