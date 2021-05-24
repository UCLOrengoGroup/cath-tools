/// \file
/// \brief The display_colourer class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_HPP

#include <iosfwd>
#include <memory>
#include <optional>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/display/display_colourer/detail/score_colour_handler.hpp"
#include "cath/display_colour/display_colour_gradient.hpp"
#include "cath/display_colour/display_colour_list.hpp"

// clang-format off
namespace cath { class alignment_free_display_colourer; }
namespace cath { class display_colour_spec; }
namespace cath { class display_spec; }
namespace cath { class viewer; }
namespace cath::align { class alignment; }
namespace cath::align { class alignment_context; }
namespace cath::file { class pdb_list; }
namespace cath::file { class strucs_context; }
namespace cath::sup { class superposition_context; }
// clang-format on

namespace cath {

	/// \brief TODOCUMENT
	class display_colourer {
	private:
		/// \brief Optional specification for post-modifying the colouring based on scores
		detail::score_colour_handler_opt the_score_colour_handler{ ::std::nullopt };

		/// \brief Pure virtual method with which each concrete display_colourer must define how to create a clone of itself
		[[nodiscard]] virtual std::unique_ptr<display_colourer> do_clone() const = 0;

		[[nodiscard]] virtual display_colour_spec do_get_colour_spec( const align::alignment_context & ) const = 0;

		[[nodiscard]] virtual std::string do_get_label() const = 0;

	  public:
		display_colourer() noexcept;
		explicit display_colourer(const detail::score_colour_handler &);
		virtual ~display_colourer() noexcept = default;

		display_colourer(const display_colourer &) = default;
		display_colourer(display_colourer &&) noexcept = default;
		display_colourer & operator=(const display_colourer &) = default;
		display_colourer & operator=(display_colourer &&) noexcept = default;

		[[nodiscard]] std::unique_ptr<display_colourer> clone() const;

		[[nodiscard]] const detail::score_colour_handler_opt &get_score_colour_handler_opt() const;
		[[nodiscard]] display_colour_spec                     get_colour_spec( const align::alignment_context & ) const;

		[[nodiscard]] std::string get_label() const;
	};

	/// \brief Default ctor
	///
	/// This is defined outside of the class because Clang 4.0 complains otherwise.
	/// Not quite sure about what that means.
	/// See Vittorio Romeo's StackOverflow question about it here:
	/// http://stackoverflow.com/questions/43819314/default-member-initializer-needed-within-definition-of-enclosing-class-outside
	inline display_colourer::display_colourer() noexcept = default;

	/// \brief NVI pass-through to the virtual do_get_label() method
	inline std::string display_colourer::get_label() const {
		return do_get_label();
	}

	bool has_score_colour_handler(const display_colourer &);
	const detail::score_colour_handler & get_score_colour_handler(const display_colourer &);

	std::unique_ptr<const display_colourer> get_display_colourer(const display_spec &,
	                                                             const display_colour_gradient & = make_default_colour_gradient() );

	display_colour_spec get_colour_spec(const display_colourer &,
	                                    const file::strucs_context &,
	                                    const align::alignment &);

	/// A named chunk of residues
	struct named_chunk {
		/// The name of the chunk
		::std::string name;

		/// The residues in the chunk
		residue_id_vec residues;
	};

	/// Type alias for a vector of named_chunk
	using named_chunk_vec = ::std::vector<named_chunk>;

	void colour_viewer_with_named_chunks( const named_chunk_vec &,
	                                      const ::std::string &,
	                                      viewer &,
	                                      ::std::ostream &,
	                                      const display_colour_list & = default_display_colour_list() );

	void colour_viewer(const display_colourer &,
	                   std::ostream &,
	                   viewer &,
	                   const align::alignment_context &);

	void colour_viewer(const alignment_free_display_colourer &,
	                   std::ostream &,
	                   viewer &,
	                   const file::pdb_list &,
	                   const str_vec &,
	                   const chop::region_vec_opt_vec &);

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_HPP
