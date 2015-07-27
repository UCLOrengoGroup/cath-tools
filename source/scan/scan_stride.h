/// \file
/// \brief The scan_stride class header

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

#ifndef SCAN_STRIDE_H_INCLUDED
#define SCAN_STRIDE_H_INCLUDED

#include <boost/operators.hpp>

#include "detail/scan_type_aliases.h"
#include "scan/detail/stride/rep_strider.h"

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		class scan_stride final : private boost::equality_comparable<scan_stride> {
		private:
			static constexpr index_type DEFAULT_STRIDE = 0;

			/// \brief TODOCUMENT
			detail::rep_strider query_from_strider { DEFAULT_STRIDE };

			/// \brief TODOCUMENT
			detail::rep_strider query_to_strider   { DEFAULT_STRIDE };

			/// \brief TODOCUMENT
			detail::rep_strider index_from_strider { DEFAULT_STRIDE };

			/// \brief TODOCUMENT
			detail::rep_strider index_to_strider   { DEFAULT_STRIDE };

		public:
			constexpr scan_stride(const index_type &arg_query_from_stride = DEFAULT_STRIDE,
			                      const index_type &arg_query_to_stride   = DEFAULT_STRIDE,
			                      const index_type &arg_index_from_stride = DEFAULT_STRIDE,
			                      const index_type &arg_index_to_stride   = DEFAULT_STRIDE);

			constexpr const detail::rep_strider & get_query_from_strider() const;
			constexpr const detail::rep_strider & get_query_to_strider() const;
			constexpr const detail::rep_strider & get_index_from_strider() const;
			constexpr const detail::rep_strider & get_index_to_strider() const;
		};

		/// \brief TODOCUMENT
		inline constexpr scan_stride::scan_stride(const index_type &arg_query_from_stride, ///< TODOCUMENT
		                                          const index_type &arg_query_to_stride,   ///< TODOCUMENT
		                                          const index_type &arg_index_from_stride, ///< TODOCUMENT
		                                          const index_type &arg_index_to_stride    ///< TODOCUMENT
		                                          ) : query_from_strider ( arg_query_from_stride ),
		                                              query_to_strider   ( arg_query_to_stride   ),
		                                              index_from_strider ( arg_index_from_stride ),
		                                              index_to_strider   ( arg_index_to_stride   ) {
		}

		/// \brief TODOCUMENT
		inline constexpr const detail::rep_strider & scan_stride::get_query_from_strider() const {
			return query_from_strider;
		}

		/// \brief TODOCUMENT
		inline constexpr const detail::rep_strider & scan_stride::get_query_to_strider() const {
			return query_to_strider;
		}

		/// \brief TODOCUMENT
		inline constexpr const detail::rep_strider & scan_stride::get_index_from_strider() const {
			return index_from_strider;
		}

		/// \brief TODOCUMENT
		inline constexpr const detail::rep_strider & scan_stride::get_index_to_strider() const {
			return index_to_strider;
		}

		/// \brief TODOCUMENT
		inline constexpr const index_type & get_query_from_stride(const scan_stride &arg_scan_stride ///< TODOCUMENT
		                                                          ) {
			return arg_scan_stride.get_query_from_strider().get_stride();
		}

		/// \brief TODOCUMENT
		inline constexpr const index_type & get_query_to_stride(const scan_stride &arg_scan_stride ///< TODOCUMENT
		                                                        ) {
			return arg_scan_stride.get_query_to_strider().get_stride();
		}

		/// \brief TODOCUMENT
		inline constexpr const index_type & get_index_from_stride(const scan_stride &arg_scan_stride ///< TODOCUMENT
		                                                          ) {
			return arg_scan_stride.get_index_from_strider().get_stride();
		}

		/// \brief TODOCUMENT
		inline constexpr const index_type & get_index_to_stride(const scan_stride &arg_scan_stride ///< TODOCUMENT
		                                                        ) {
			return arg_scan_stride.get_index_to_strider().get_stride();
		}

		/// \brief TODOCUMENT
		inline constexpr bool operator==(const scan_stride &arg_scan_stride_a, ///< TODOCUMENT
		                                 const scan_stride &arg_scan_stride_b  ///< TODOCUMENT
		                                 ) {
			return (
				arg_scan_stride_a.get_query_from_strider() == arg_scan_stride_b.get_query_from_strider()
				&&
				arg_scan_stride_a.get_query_to_strider()   == arg_scan_stride_b.get_query_to_strider()
				&&
				arg_scan_stride_a.get_index_from_strider() == arg_scan_stride_b.get_index_from_strider()
				&&
				arg_scan_stride_a.get_index_to_strider()   == arg_scan_stride_b.get_index_to_strider()
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates scan_stride
		inline constexpr index_type from_co_stride(const scan_stride &arg_stride ///< The scan_stride object containing the relevant stride values
		                                           ) {
			return co_stride(
				arg_stride.get_query_from_strider(),
				arg_stride.get_index_from_strider()
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates scan_stride
		inline constexpr index_type to_co_stride(const scan_stride &arg_stride ///< The scan_stride object containing the relevant stride values
		                                         ) {
			return co_stride(
				arg_stride.get_query_to_strider(),
				arg_stride.get_index_to_strider()
			);
		}

		namespace detail {
			rep_rep_pair_opt get_from_rep_of_indices(const scan_stride &,
			                                         const index_type &,
			                                         const index_type &);
			rep_rep_pair_opt get_to_rep_of_indices(const scan_stride &,
			                                       const index_type &,
			                                       const index_type &);
		}

	}
}

#endif
