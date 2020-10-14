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

#ifndef _CATH_TOOLS_SOURCE_UNI_DISPLAY_VIEWER_VIEWER_HPP
#define _CATH_TOOLS_SOURCE_UNI_DISPLAY_VIEWER_VIEWER_HPP

#include <string>

#include <boost/utility/string_ref.hpp>
#include <boost/filesystem/path.hpp>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/display/colour_category.hpp"

namespace cath { class display_colour; }
namespace cath { class display_colourer; }
namespace cath { class display_spec; }
namespace cath { namespace align { class alignment; } }
namespace cath { namespace align { class alignment_context; } }
namespace cath { namespace file { class pdb_list; } }
namespace cath { namespace sup { class superposition; } }
namespace cath { namespace sup { class superposition_content_spec; } }
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
		virtual bool do_accepts_multiple_colourings() const;
		virtual void do_begin_colouring(std::ostream &);
		virtual std::string do_get_colour_base_str(const std::string &) const = 0;
		virtual std::string do_get_colour_pdb_str(const std::string &,
		                                          const std::string &) const = 0;
		virtual std::string do_get_colour_pdb_residues_str(const std::string &,
		                                                   const std::string &,
		                                                   const residue_id_vec &) const = 0;
		virtual void do_end_colouring(::std::ostream &,
		                              const ::std::string &);
		virtual void do_write_alignment_extras(std::ostream &,
		                                       const sup::superposition_context &) const = 0;
		virtual void do_write_end(std::ostream &,
		                          const boost::string_ref &) const = 0;

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
		bool accepts_multiple_colourings() const;
		void begin_colouring(std::ostream &);
		std::string get_colour_base_str(const std::string &) const;
		std::string get_colour_pdb_str(const std::string &,
		                               const std::string &) const;
		std::string get_colour_pdb_residues_str(const std::string &,
		                                        const std::string &,
		                                        const residue_id_vec &) const;
		void end_colouring(std::ostream &,
		                   const ::std::string &);
		void write_alignment_extras(std::ostream &,
		                            const sup::superposition_context &) const;
		void write_end(std::ostream &,
		               const boost::string_ref &) const;
	};

	std::string clean_name_for_viewer(const std::string &);

	str_vec clean_names_for_viewer(const str_vec &);
	
	str_vec clean_names_for_viewer(const sup::superposition_context &);

	str_vec clean_names_for_viewer(const align::alignment_context &);

	/// \brief Policy for handling a missing alignment where a display_spec requires one for the colouring
	enum class missing_aln_policy : bool {
		WARN_AND_COLOUR_CONSECUTIVELY, ///< Warn on encountering a missing alignment and just use consecutive colouring
		THROW                          ///< Throw on encountering a missing alignment
	};

	void output_superposition_to_viewer( std::ostream &,
	                                     viewer &,
	                                     const display_spec &,
	                                     const sup::superposition_context &,
	                                     const sup::superposition_content_spec &,
	                                     const missing_aln_policy & = missing_aln_policy::THROW );

	void output_superposition_to_viewer_file( const boost::filesystem::path &,
	                                          viewer &,
	                                          const display_spec &,
	                                          const sup::superposition_context &,
	                                          const sup::superposition_content_spec &,
	                                          const missing_aln_policy & = missing_aln_policy::THROW );

	std::string colour_of_index_from_colours_string(const size_t &,
	                                                const std::string &);

	std::string base_colour_name();

	std::string generate_colour_name(const size_t &,
	                                 const size_t &,
	                                 const colour_category &);

	str_vec generate_colour_names(const size_t &,
	                              const colour_category &);

	/// \brief NVI pass-through to the virtual do_accepts_multiple_colourings()
	inline bool viewer::accepts_multiple_colourings() const {
		return do_accepts_multiple_colourings();
	}

	/// \brief NVI pass-through to the virtual do_begin_colouring()
	inline void viewer::begin_colouring(std::ostream &prm_os ///< The ostream to which the PyMOL commands should be written
	                                    ) {
		do_begin_colouring( prm_os );
	}

	/// \brief NVI pass-through to the virtual do_end_colouring()
	///
	/// \param prm_os    The ostream to which the PyMOL commands should be written
	/// \param prm_label The label to be used for the colouring that is ending
	inline void viewer::end_colouring( ::std::ostream &prm_os, const ::std::string &prm_label ) {
		do_end_colouring( prm_os, prm_label );
	}

} // namespace cath

#endif
