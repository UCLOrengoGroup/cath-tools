/// \file
/// \brief The score type_aliases header

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

#ifndef SCORE_TYPE_ALIASES_H_INCLUDED
#define SCORE_TYPE_ALIASES_H_INCLUDED

#include <boost/config.hpp> /// \todo Come a resolution for Boost Trac tickets 12142 & 12179, remove this #include
#include <boost/optional/optional_fwd.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <cstddef>
#include <vector>

namespace cath { namespace homcheck { class ssap_and_prc; } }
namespace cath { namespace score { class aligned_pair_score; } }
namespace cath { namespace score { class classn_stat_pair_series; } }
namespace cath { namespace score { class named_true_false_pos_neg_list; } }
namespace cath { namespace score { class score_classn_value; } }
namespace cath { namespace score { class score_classn_value_list; } }
namespace cath { namespace score { class substitution_matrix; } }
namespace cath { namespace score { class true_false_pos_neg; } }
namespace cath { namespace score { class value_list_scaling; } }

namespace cath {
	namespace homcheck {
		/// \brief Type alias for a vector of ssap_and_prc objects
		using ssap_and_prc_vec      = std::vector<ssap_and_prc>;

		/// \brief Type alias for a reference_wrapper to a const ssap_and_prc
		using ssap_and_prc_cref     = std::reference_wrapper<const ssap_and_prc>;

		/// \brief Type alias for an optional reference_wrapper to a const ssap_and_prc
		using ssap_and_prc_cref_opt = boost::optional<ssap_and_prc_cref>;
	}
}

namespace cath {
	namespace score {
		/// \brief TODOCUMENT
		using str_aligned_pair_score_pmap = boost::ptr_map<std::string, aligned_pair_score>;


		/// \brief TODOCUMENT
		using classn_stat_pair_series_vec      = std::vector<classn_stat_pair_series>;

		/// \brief TODOCUMENT
		using classn_stat_pair_series_vec_itr  = classn_stat_pair_series_vec::iterator;

		/// \brief TODOCUMENT
		using classn_stat_pair_series_vec_citr = classn_stat_pair_series_vec::const_iterator;


		/// \brief TODOCUMENT
		using named_true_false_pos_neg_list_vec      = std::vector<named_true_false_pos_neg_list>;

		/// \brief TODOCUMENT
		using named_true_false_pos_neg_list_vec_itr  = named_true_false_pos_neg_list_vec::iterator;

		/// \brief TODOCUMENT
		using named_true_false_pos_neg_list_vec_citr = named_true_false_pos_neg_list_vec::const_iterator;


		/// \brief TODOCUMENT
		using substitution_matrix_vec = std::vector<substitution_matrix>;


		/// \brief TODOCUMENT
		using true_false_pos_neg_vec      = std::vector<true_false_pos_neg>;

		/// \brief TODOCUMENT
		using true_false_pos_neg_vec_citr = true_false_pos_neg_vec::const_iterator;


		/// \brief TODOCUMENT
		using doub_true_false_pos_neg_pair = std::pair<double, true_false_pos_neg>;

		/// \brief TODOCUMENT
		using doub_true_false_pos_neg_pair_vec = std::vector<doub_true_false_pos_neg_pair>;


		/// \brief TODOCUMENT
		using score_classn_value_vec      = std::vector<score_classn_value>;

		/// \brief TODOCUMENT
		using score_classn_value_vec_itr  = score_classn_value_vec::iterator;

		/// \brief TODOCUMENT
		using score_classn_value_vec_citr = score_classn_value_vec::const_iterator;


		/// \brief TODOCUMENT
		using score_classn_value_vec_vec = std::vector<score_classn_value_vec>;


		/// \brief TODOCUMENT
		using score_classn_value_list_vec      = std::vector<score_classn_value_list>;

		/// \brief TODOCUMENT
		using score_classn_value_list_vec_itr  = score_classn_value_list_vec::iterator;

		/// \brief TODOCUMENT
		using score_classn_value_list_vec_citr = score_classn_value_list_vec::const_iterator;


		/// \brief The type to be used for the calculated scores
		using score_value = double;

		/// \brief A convenience type-alias for a vector of score_values
		using score_value_vec = std::vector<score_value>;

		/// \brief TODOCUMENT
		using score_size_pair = std::pair<score_value, size_t>;

//		class sensitivity;
//		class specificity;
//		class precision;
//		class fall_out;

//		/// \brief A type alias permitting sensitivity to also be called recall
//		using recall = sensitivity;
//
//		/// \brief A type alias permitting sensitivity to also be called hit_rate
//		using hit_rate = sensitivity;
//
//		/// \brief A type alias permitting sensitivity to also be called true_positive_rate
//		using true_positive_rate = sensitivity;
//
//		/// \brief A type alias permitting specificity to also be called true_negative_rate
//		using true_negative_rate = specificity;
//
//		/// \brief A type alias permitting precision to also be called positive_predictive_value
//		using positive_predictive_value = precision;
//
//		/// \brief A type alias permitting fall_out to also be called false_positive_rate
//		using false_positive_rate = fall_out;

		namespace detail {
			class score_common_coord_handler;

			/// \brief TODOCUMENT
			using score_common_coord_handler_vec = std::vector<score_common_coord_handler>;
		}

		/// \brief TODOCUMENT
		using value_list_scaling_vec = std::vector<value_list_scaling>;
	}
}
#endif
