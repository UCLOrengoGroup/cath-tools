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
			superposition_context(const superposition &,
			                      const file::strucs_context &);

			superposition_context(const superposition &,
			                      const file::strucs_context &,
			                      const align::alignment &);

			superposition_context(const superposition &,
			                      const file::pdb_list &,
			                      const str_vec &,
			                      const chop::region_vec_opt_vec &);

			superposition_context(const superposition &,
			                      const file::pdb_list &,
			                      const str_vec &,
			                      const chop::region_vec_opt_vec &,
			                      const align::alignment &);

			const superposition & get_superposition() const;
			const file::strucs_context & get_strucs_context() const;
			bool has_alignment() const;
			const align::alignment & get_alignment() const;

			superposition_context & set_pdbs(const file::pdb_list &);
		};

		const file::pdb_list & get_pdbs(const superposition_context &);
		const str_vec & get_names(const superposition_context &);
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

		align::alignment_context make_alignment_context(const superposition_context &);

		superposition_context superposition_context_from_ptree(const boost::property_tree::ptree &);

		void save_to_ptree(boost::property_tree::ptree &,
		                   const superposition_context &);

		boost::property_tree::ptree make_ptree_of(const superposition_context &);

		superposition_context superposition_context_from_json_string(const std::string &);

		std::string to_json_string(const superposition_context &,
		                           const common::json_style & = common::json_style::PRETTY);

		superposition_context read_superposition_context_from_json_file(const boost::filesystem::path &);

		void write_to_json_file(const boost::filesystem::path &,
		                        const superposition_context &,
		                        const common::json_style & = common::json_style::PRETTY);
	} // namespace sup
} // namespace cath

#endif
