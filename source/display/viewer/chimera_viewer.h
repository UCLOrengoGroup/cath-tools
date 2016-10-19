/// \file
/// \brief The chimera_viewer class header

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

#ifndef _CATH_TOOLS_SOURCE_DISPLAY_VIEWER_CHIMERA_VIEWER_H
#define _CATH_TOOLS_SOURCE_DISPLAY_VIEWER_CHIMERA_VIEWER_H

#include "display/viewer/viewer.h"

namespace cath {
	namespace detail {
		void write_chimera_pair_alignments(std::ostream &,
		                                   const sup::superposition_context &);

		void write_chimera_global_alignment(std::ostream &,
		                                    const sup::superposition_context &);
	} // namespace detail

	/// \brief TODOCUMENT
	class chimera_viewer final : public viewer {
	private:
		virtual std::string do_default_executable() const override final;
		virtual std::string do_default_file_extension() const override final;

		virtual void do_write_start(std::ostream &) const override final;
		virtual void do_write_load_pdbs(std::ostream &,
		                                const sup::superposition &,
		                                const file::pdb_list &,
		                                const str_vec &) const override final;
		virtual void do_define_colour(std::ostream &,
		                              const display_colour &,
		                              const std::string &) const override final;
		virtual void do_colour_pdb(std::ostream &,
		                           const std::string &,
		                           const std::string &) const override final;
		virtual void do_colour_pdb_residues(std::ostream &,
		                                    const std::string &,
		                                    const std::string &,
		                                    const residue_name_vec &) const override final;
		virtual void do_write_alignment_extras(std::ostream &,
		                                       const sup::superposition_context &) const override final;
		virtual void do_write_end(std::ostream &) const override final;

	public:
		virtual ~chimera_viewer() noexcept = default;

		static std::string parse_residue_name_for_chimera(const residue_name &);
		static str_vec parse_residue_names_for_chimera(const residue_name_vec &);

		/// \brief TODOCUMENT
		static constexpr size_t RESIDUE_BATCH_SIZE = 200;
	};



} // namespace cath

#endif
