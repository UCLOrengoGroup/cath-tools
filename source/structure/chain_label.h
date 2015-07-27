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

#ifndef CHAIN_LABEL_H_INCLUDED
#define CHAIN_LABEL_H_INCLUDED

#include <boost/operators.hpp>

#include <iosfwd>

namespace cath {

	/// \brief TODOCUMENT
	///
	/// Can much of this be made constexpr?
	class chain_label final : private boost::equality_comparable<chain_label> {
	private:
		friend bool operator==(const chain_label &, const chain_label &);

		/// \brief TODOCUMENT
		char chain_char;

		void sanity_check() const;

		const char & get_char() const;

	public:
		explicit chain_label(const char &);
		
		std::string to_string() const;
	};

	/// \brief TODOCUMENT
	inline const char & chain_label::get_char() const {
		return chain_char;
	}

	bool operator==(const chain_label &,
	                const chain_label &);

	std::ostream & operator<<(std::ostream &,
	                          const chain_label &);
}

#endif
