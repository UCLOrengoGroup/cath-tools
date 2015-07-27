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

#ifndef DOMAIN_H_INCLUDED
#define DOMAIN_H_INCLUDED

#include <boost/optional.hpp>

#include "chopping/chopping_type_aliases.h"
#include "chopping/region/region.h"
#include "common/type_aliases.h"

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		///
		/// Invariants:
		///  * all segments must have the same residue_locating
		///     (ie whether they locate their residues by names and/or indices)
		class domain final {
		private:
			/// \brief TODOCUMENT
			region_vec segments;

			/// \brief TODOCUMENT
			opt_str domain_id;

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

			void set_opt_domain_id(const opt_str &);
			const opt_str & get_opt_domain_id() const;

//			iterator begin();
//			iterator end();
			const_iterator begin() const;
			const_iterator end() const;
		};

		bool has_domain_id(const domain &);
		std::string get_domain_id(const domain &);

		opt_residue_locating get_residue_locating(const domain &);

	}
}

#endif
