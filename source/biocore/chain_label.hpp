/// \file
/// \brief The chain_label class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_CHAIN_LABEL_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_CHAIN_LABEL_H

#include <boost/operators.hpp>

#include <iosfwd>

namespace cath {

	/// \brief TODOCUMENT
	class chain_label final : private boost::equality_comparable<chain_label> {
	private:
		friend constexpr bool operator==(const chain_label &,
		                                 const chain_label &);
		friend constexpr bool operator<(const chain_label &,
		                                const chain_label &);

		/// \brief TODOCUMENT
		char chain_char = 0;

		constexpr const char & get_char() const;

	public:
		constexpr chain_label() noexcept = default;
		explicit constexpr chain_label(const char &);

		bool is_null() const;

		std::string to_string() const;
	};

	/// \brief TODOCUMENT
	inline constexpr const char & chain_label::get_char() const {
		return chain_char;
	}

	/// \brief Ctor for chain_label
	inline constexpr chain_label::chain_label(const char &arg_chain_char ///< TODOCUMENT
	                                          ) : chain_char( arg_chain_char ) {
	}

	/// \brief TODOCUMENT
	inline bool chain_label::is_null() const {
		return ( get_char() == 0 );
	}

	/// \brief TODOCUMENT
	inline constexpr bool operator==(const chain_label &arg_chain_label_a, ///< TODOCUMENT
	                                 const chain_label &arg_chain_label_b  ///< TODOCUMENT
	                                 ) {
		return ( arg_chain_label_a.get_char() == arg_chain_label_b.get_char() );
	}

	/// \brief TODOCUMENT
	inline constexpr bool operator<(const chain_label &arg_chain_label_a, ///< TODOCUMENT
	                                const chain_label &arg_chain_label_b  ///< TODOCUMENT
	                                ) {
		return ( arg_chain_label_a.get_char() < arg_chain_label_b.get_char() );
	}

	std::ostream & operator<<(std::ostream &,
	                          const chain_label &);
} // namespace cath

#endif
