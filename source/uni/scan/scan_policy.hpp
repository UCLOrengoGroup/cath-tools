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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCAN_SCAN_POLICY_H
#define _CATH_TOOLS_SOURCE_UNI_SCAN_SCAN_POLICY_H

#include "scan/quad_criteria.hpp"
#include "scan/res_pair_keyer/res_pair_keyer.hpp"
#include "scan/scan_stride.hpp"

namespace cath {
	namespace scan {

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
			const quad_criteria & get_criteria() const;
			const scan_stride & get_scan_stride() const;
		};

		/// \brief TODOCUMENT
		template <typename... KPs>
		scan_policy<KPs...>::scan_policy(const res_pair_keyer<KPs...> &arg_keyer,      ///< TODOCUMENT
		                                 quad_criteria                 arg_criteria,   ///< TODOCUMENT
		                                 scan_stride                   arg_scan_stride ///< TODOCUMENT
		                                 ) : keyer        { arg_keyer                    },
		                                     the_criteria { std::move( arg_criteria    ) },
		                                     the_stride   { std::move( arg_scan_stride ) } {
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
		scan_policy<KPs...> make_scan_policy(const res_pair_keyer<KPs...> &arg_keyer,      ///< TODOCUMENT
		                                     const quad_criteria       &arg_criteria,   ///< TODOCUMENT
		                                     const scan_stride         &arg_scan_stride ///< TODOCUMENT
		                                     ) {
			return {
				arg_keyer,
				arg_criteria,
				arg_scan_stride
			};
		}

	} // namespace scan
} // namespace cath

#endif
