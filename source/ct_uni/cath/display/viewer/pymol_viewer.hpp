/// \file
/// \brief The pymol_viewer class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_VIEWER_PYMOL_VIEWER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_VIEWER_PYMOL_VIEWER_HPP

#include "cath/display/viewer/viewer.hpp"

namespace cath {
	namespace detail {
		void write_pymol_pair_alignments(std::ostream &,
		                                 const sup::superposition_context &);

		void write_pymol_global_alignment(std::ostream &,
		                                  const sup::superposition_context &);
	} // namespace detail

	/// \brief TODOCUMENT
	class pymol_viewer final : public viewer {
	private:
		/// \brief A counter for the scenes that get stored
		///
		/// This implies the F key (eg F1, F2, ...) under which the scene is stored
		size_t scene_count = 1;

		void record_scene(std::ostream &,
		                  const std::string &);

		std::string do_default_executable() const final;
		std::string do_default_file_extension() const final;

		void do_write_start(std::ostream &) const final;
		void do_write_load_pdbs(std::ostream &,
		                        const sup::superposition &,
		                        const file::pdb_list &,
		                        const str_vec &) const final;
		void do_define_colour(std::ostream &,
		                      const display_colour &,
		                      const std::string &) const final;
		bool do_accepts_multiple_colourings() const final;

		void do_begin_colouring(std::ostream &) final;
		std::string do_get_colour_base_str(const std::string &) const final;
		std::string do_get_colour_pdb_str(const std::string &,
		                                  const std::string &) const final;
		std::string do_get_colour_pdb_residues_str(const std::string &,
		                                           const std::string &,
		                                           const residue_id_vec &) const final;
		void do_end_colouring(std::ostream &,
		                      const ::std::string &) final;

		void do_write_alignment_extras(std::ostream &,
		                               const sup::superposition_context &) const final;
		void do_write_end(std::ostream &,
		                  const boost::string_ref &) const final;

	public:
		/// \brief TODOCUMENT
		static constexpr size_t RESIDUE_BATCH_SIZE = 200;
	};



} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_VIEWER_PYMOL_VIEWER_HPP
