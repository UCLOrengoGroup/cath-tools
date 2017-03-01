/// \file
/// \brief The domain class header

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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_DOMAIN_DOMAIN_H
#define _CATH_TOOLS_SOURCE_CHOPPING_DOMAIN_DOMAIN_H

#include <boost/operators.hpp>
#include <boost/optional.hpp>

#include "chopping/chopping_type_aliases.hpp"
#include "chopping/region/region.hpp"
#include "common/type_aliases.hpp"

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		///
		/// Invariants:
		///  * all segments must have the same residue_locating
		///     (ie whether they locate their residues by names and/or indices)
		class domain final : private boost::equality_comparable<domain>  {
		private:
			/// \brief TODOCUMENT
			region_vec segments;

			/// \brief TODOCUMENT
			str_opt domain_id;

			void sanity_check() const;

		public:
//			using iterator = region_vec::iterator;
			using const_iterator = region_vec::const_iterator;

			explicit domain(const region_vec &);
			explicit domain(const region_vec &,
			                const std::string &);

			size_t num_segments() const;
//			region operator[](const size_t &);
			const region & operator[](const size_t &) const;

			void set_opt_domain_id(const str_opt &);
			const str_opt & get_opt_domain_id() const;

//			iterator begin();
//			iterator end();
			const_iterator begin() const;
			const_iterator end() const;
		};

		bool operator==(const domain &,
		                const domain &);

		bool has_domain_id(const domain &);
		std::string get_domain_id(const domain &);

		residue_locating_opt get_residue_locating(const domain &);

		std::string to_string(const domain &);
		std::ostream & operator<<(std::ostream &,
		                          const domain &);

	} // namespace chop
} // namespace cath

#endif
