/// \file
/// \brief The value_list_scaling class header

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

#ifndef VALUE_LIST_SCALING_H_INCLUDED
#define VALUE_LIST_SCALING_H_INCLUDED

#include <iosfwd>

namespace cath {
	namespace score {

		/// \brief Represent a particular linear scaling of the form f(x) = mx + c
		class value_list_scaling final {
		private:
			/// \brief TODOCUMENT
			double multiplier = 1.0;

			/// \brief TODOCUMENT
			double constant = 0.0;

		public:
			value_list_scaling(const double &,
			                   const double &);

			const double & get_multiplier() const;
			const double & get_constant() const;
		};

		std::string to_string(const value_list_scaling &);

		std::ostream & operator<<(std::ostream &,
		                          const value_list_scaling &);

		void scale_value(const value_list_scaling &,
		                 double &);

		double scale_value_copy(const value_list_scaling &,
		                        double);

	}
}

#endif

