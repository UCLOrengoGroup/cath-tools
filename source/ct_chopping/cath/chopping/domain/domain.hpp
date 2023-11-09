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

#ifndef CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_DOMAIN_DOMAIN_HPP
#define CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_DOMAIN_DOMAIN_HPP

#include <cstddef>
#include <ostream>
#include <string>

#include <boost/operators.hpp>

#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath::chop {

	/// \brief TODOCUMENT
	///
	/// Invariants:
	///  * all segments must have the same residue_locating
	///     (ie whether they locate their residues by names and/or indices)
	///
	/// \todo Should a domain be allowed to have no segments and if so, should that mean
	///       that the domain covers everything (eg the whole PDB)?
	class domain final : private boost::equality_comparable<domain>  {
	private:
		/// \brief TODOCUMENT
		region_vec segments;

		/// \brief TODOCUMENT
		str_opt domain_id;

		void sanity_check() const;

	public:
//		using iterator = region_vec::iterator;
		using const_iterator = region_vec::const_iterator;

		/// \brief Type alias for the value_type
		///
		/// This is added because otherwise the Travis-CI Mac build fails to compile
		/// due to a use of value_type in `boost/test/utils/is_forward_iterable.hpp`
		/// that gets invoked by the `BOOST_TEST( [...] == [...] )` comparisons of
		/// domain in sillitoe_chopping_format_test.cpp
		using value_type     = region_vec::value_type;

		explicit domain(region_vec);
		explicit domain(region_vec,
		                std::string);

		[[nodiscard]] size_t num_segments() const;
		// region operator[](const size_t &);
		const region & operator[](const size_t &) const;

		void set_opt_domain_id(const str_opt &);
		[[nodiscard]] const str_opt &get_opt_domain_id() const;

		// iterator begin();
		// iterator end();
		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	region_vec_opt get_regions_opt(const domain_opt &);

	bool operator==(const domain &,
	                const domain &);

	bool has_domain_id(const domain &);
	std::string get_domain_id(const domain &);

	residue_locating_opt get_residue_locating(const domain &);

	std::string to_string(const domain &);
	std::ostream & operator<<(std::ostream &,
	                          const domain &);

} // namespace cath::chop

#endif // CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_DOMAIN_DOMAIN_HPP
