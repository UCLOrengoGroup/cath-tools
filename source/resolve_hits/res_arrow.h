/// \file
/// \brief The res_arrow class header

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

#ifndef RES_ARROW_H_INCLUDED
#define RES_ARROW_H_INCLUDED

#include <boost/operators.hpp>

#include "resolve_hits/resolve_hits_type_aliases.h"

#include <iosfwd>

namespace cath {
	namespace rslv {

		/// \brief TODOCUMENT
		class res_arrow final : private boost::totally_ordered<res_arrow,
		                                boost::additive<res_arrow, resarw_t> > {
		private:
			friend constexpr res_arrow arrow_before_res(const residx_t &);
			friend constexpr res_arrow arrow_after_res (const residx_t &);

			/// \brief TODOCUMENT
			resarw_t arrow;

			explicit constexpr res_arrow(const resarw_t &);

		public:
			constexpr residx_t res_before() const;
			constexpr const residx_t & res_after () const;

			constexpr const resarw_t & get_index() const;

			res_arrow & operator+=(const resarw_t &);
			res_arrow & operator-=(const resarw_t &);
		};

		constexpr bool operator<(const res_arrow &arg_arrow_a, ///< TODOCUMENT
		                         const res_arrow &arg_arrow_b  ///< TODOCUMENT
		                         ) {
			return arg_arrow_a.get_index() < arg_arrow_b.get_index();
		}

		constexpr bool operator==(const res_arrow &arg_arrow_a, ///< TODOCUMENT
		                          const res_arrow &arg_arrow_b  ///< TODOCUMENT
		                          ) {
			return arg_arrow_a.get_index() == arg_arrow_b.get_index();
		}

		std::string to_string(const res_arrow &);
		std::ostream & operator<<(std::ostream &,
		                          const res_arrow &);

		constexpr res_arrow arrow_before_res(const residx_t &);
		constexpr res_arrow arrow_after_res (const residx_t &);
		constexpr res_arrow start_arrow();

		/// \brief TODOCUMENT
		inline constexpr res_arrow::res_arrow(const resarw_t &arg_res_arrow_index ///< TODOCUMENT
		                                      ) : arrow ( arg_res_arrow_index ) {
		}

		/// \brief TODOCUMENT
		inline constexpr residx_t res_arrow::res_before() const {
#ifndef NDEBUG
			return ( arrow > 0 ) ? ( arrow - 1 ) : throw( "Cannot return a residue before an arrow that precedes the first residue" );
#else
			return                 ( arrow - 1 );
#endif
		}

		/// \brief TODOCUMENT
		inline constexpr const residx_t & res_arrow::res_after () const {
			return arrow;
		}

		/// \brief TODOCUMENT
		inline constexpr const resarw_t & res_arrow::get_index() const {
			return arrow;
		}

		/// \brief TODOCUMENT
		inline res_arrow & res_arrow::operator+=(const resarw_t &arg_offset ///< TODOCUMENT
		                                         ) {
			arrow += arg_offset;
			return *this;
		}

		/// \brief TODOCUMENT
		inline res_arrow & res_arrow::operator-=(const resarw_t &arg_offset ///< TODOCUMENT
		                                         ) {
			arrow -= arg_offset;
			return *this;
		}

		/// \brief TODOCUMENT
		inline constexpr res_arrow arrow_before_res(const residx_t &arg_res_index ///< TODOCUMENT
		                                            ) {
			/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
			return res_arrow{ arg_res_index     };
		}

		/// \brief TODOCUMENT
		inline constexpr res_arrow arrow_after_res(const residx_t &arg_res_index ///< TODOCUMENT
		                                           ) {
			/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
			return res_arrow{ arg_res_index + 1 };
		}

		/// \brief TODOCUMENT
		inline constexpr res_arrow start_arrow() {
			return arrow_before_res( 0 );
		}



	}
}

#endif
