/// \file
/// \brief The scan_policy class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SCAN_POLICY_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SCAN_POLICY_HPP

#include "cath/scan/quad_criteria.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer.hpp"
#include "cath/scan/scan_stride.hpp"

namespace cath::scan {

	/// \brief TODOCUMENT
	///
	/// It is helpful that this is read-only after construction because that makes is certain that
	/// anything using the scan_policy will have been using consistent values over time.
	///
	/// It also ensures that  a scan_index and scan_query_set using the same scan_policy object must
	/// have been using the same values as each other.
	///
	/// \todo Move the strides into a subclass that can be passed around independently without the need for a template parameter
	template <typename... KPs>
	class scan_policy final {
	private:
		/// \brief TODOCUMENT
		res_pair_keyer<KPs...> keyer;

		/// \brief TODOCUMENT
		quad_criteria the_criteria;

		/// \brief TODOCUMENT
		scan_stride the_stride;

	public:
		scan_policy(const res_pair_keyer<KPs...> &,
		            quad_criteria,
		            scan_stride);

		const res_pair_keyer<KPs...> & get_keyer() const;
		[[nodiscard]] const quad_criteria &get_criteria() const;
		[[nodiscard]] const scan_stride &  get_scan_stride() const;
	};

	/// \brief TODOCUMENT
	template <typename... KPs>
	scan_policy<KPs...>::scan_policy(const res_pair_keyer<KPs...> &prm_keyer,      ///< TODOCUMENT
	                                 quad_criteria                 prm_criteria,   ///< TODOCUMENT
	                                 scan_stride                   prm_scan_stride ///< TODOCUMENT
	                                 ) : keyer        { prm_keyer                    },
	                                     the_criteria { std::move( prm_criteria    ) },
	                                     the_stride   { std::move( prm_scan_stride ) } {
	}

	/// \brief TODOCUMENT
	template <typename... KPs>
	const res_pair_keyer<KPs...> & scan_policy<KPs...>::get_keyer() const {
		return keyer;
	}

	/// \brief TODOCUMENT
	template <typename... KPs>
	const quad_criteria & scan_policy<KPs...>::get_criteria() const {
			return the_criteria;
	}

	/// \brief TODOCUMENT
	template <typename... KPs>
	const scan_stride & scan_policy<KPs...>::get_scan_stride() const {
		return the_stride;
	}

	/// \brief TODOCUMENT
	template <typename... KPs>
	scan_policy<KPs...> make_scan_policy(const res_pair_keyer<KPs...> &prm_keyer,      ///< TODOCUMENT
	                                     const quad_criteria       &prm_criteria,   ///< TODOCUMENT
	                                     const scan_stride         &prm_scan_stride ///< TODOCUMENT
	                                     ) {
		return {
			prm_keyer,
			prm_criteria,
			prm_scan_stride
		};
	}

} // namespace cath::scan

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SCAN_POLICY_HPP
