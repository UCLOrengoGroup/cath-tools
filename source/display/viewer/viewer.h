/// \file
/// \brief The viewer class header

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

#ifndef _CATH_TOOLS_SOURCE_DISPLAY_VIEWER_VIEWER_H
#define _CATH_TOOLS_SOURCE_DISPLAY_VIEWER_VIEWER_H

#include "common/type_aliases.h"
#include "structure/structure_type_aliases.h"

#include <string>

namespace cath { namespace align { class alignment; } }
namespace cath { namespace align { class alignment_context; } }
namespace cath { class display_colour; }
namespace cath { class display_spec; }
namespace cath { namespace file { class pdb_list; } }
namespace cath { namespace sup { class superposition; } }
namespace cath { namespace sup { class superposition_context; } }

namespace cath {
	/// \brief TODOCUMENT
	class viewer {
	private:
		virtual std::string do_default_executable() const = 0;
		virtual std::string do_default_file_extension() const = 0;

		virtual void do_write_start(std::ostream &) const = 0;
		virtual void do_write_load_pdbs(std::ostream &,
		                                const sup::superposition &,
		                                const file::pdb_list &,
		                                const str_vec &) const = 0;
		virtual void do_define_colour(std::ostream &,
		                              const display_colour &,
		                              const std::string &) const = 0;
		virtual std::string do_get_colour_pdb_str(const std::string &,
		                                          const std::string &) const = 0;
		virtual std::string do_get_colour_pdb_residues_str(const std::string &,
		                                                   const std::string &,
		                                                   const residue_name_vec &) const = 0;
		virtual void do_write_alignment_extras(std::ostream &,
		                                       const sup::superposition_context &) const = 0;
		virtual void do_write_end(std::ostream &) const = 0;

	public:
		viewer() = default;
		virtual ~viewer() noexcept = default;

		viewer(const viewer &) = default;
		viewer(viewer &&) noexcept = default;
		viewer & operator=(const viewer &) = default;
		viewer & operator=(viewer &&) noexcept = default;

		std::string default_executable() const;
		std::string default_file_extension() const;

		void write_start(std::ostream &) const;
		void write_load_pdbs(std::ostream &,
		                     const sup::superposition &,
		                     const file::pdb_list &,
		                     const str_vec &) const;
		void define_colour(std::ostream &,
		                   const display_colour &,
		                   const std::string &) const;
		std::string get_colour_pdb_str(const std::string &,
		                               const std::string &) const;
		std::string get_colour_pdb_residues_str(const std::string &,
		                                        const std::string &,
		                                        const residue_name_vec &) const;
		void write_alignment_extras(std::ostream &,
		                            const sup::superposition_context &) const;
		void write_end(std::ostream &) const;
	};

	std::string clean_name_for_viewer(const std::string &);

	str_vec clean_names_for_viewer(const str_vec &);
	
	str_vec clean_names_for_viewer(const sup::superposition_context &);

	str_vec clean_names_for_viewer(const align::alignment_context &);

	void output_superposition_to_viewer(std::ostream &,
	                                    const viewer &,
	                                    const display_spec &,
	                                    const sup::superposition_context &,
	                                    const bool & = false);

	std::string colour_of_index_from_colours_string(const size_t &,
	                                                const std::string &);

	std::string generate_colour_name(const size_t &,
	                                 const size_t &);

	str_vec generate_colour_names(const size_t &);

} // namespace cath

#endif
