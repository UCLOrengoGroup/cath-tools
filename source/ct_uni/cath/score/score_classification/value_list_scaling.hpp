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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_VALUE_LIST_SCALING_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_VALUE_LIST_SCALING_HPP

#include <iosfwd>

namespace cath::score {

	/// \brief Represent a particular linear scaling of the form f(x) = mx + c
	class value_list_scaling final {
	private:
		/// \brief TODOCUMENT
		double multiplier = 1.0;

		/// \brief TODOCUMENT
		double constant = 0.0;

	public:
		constexpr value_list_scaling(const double &,
		                             const double &);

		[[nodiscard]] constexpr const double &get_multiplier() const;
		[[nodiscard]] constexpr const double &get_constant() const;

		/// \brief The value to assign to bad (probably absent) entries after scaling
		static constexpr double BAD_SCALED_VALUE = -999.0;
	};

	/// \brief Ctor from the multiplier (m) and constant (c)
	inline constexpr value_list_scaling::value_list_scaling(const double &prm_multiplier, ///< The multiplier by which values should be multiplied
	                                                        const double &prm_constant    ///< The constant by which the multiplied values should be increased
	                                                        ) : multiplier ( prm_multiplier ),
	                                                            constant   ( prm_constant   ) {
	}

	/// \brief Getter for the multiplier
	inline constexpr const double & value_list_scaling::get_multiplier() const {
		return multiplier;
	}

	/// \brief Getter for the constant
	inline constexpr const double & value_list_scaling::get_constant() const {
		return constant;
	}

	std::string to_string(const value_list_scaling &);

	std::ostream & operator<<(std::ostream &,
	                          const value_list_scaling &);

	void scale_value(const value_list_scaling &,
	                 double &);

	double scale_value_copy(const value_list_scaling &,
	                        double);

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_VALUE_LIST_SCALING_HPP

