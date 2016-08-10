/// \file
/// \brief The score_value_reader class header

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

#ifndef SCORE_VALUE_READER_H_INCLUDED
#define SCORE_VALUE_READER_H_INCLUDED

#include <boost/filesystem/path.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

namespace cath { namespace score { class aligned_pair_score_value_list; } }

namespace cath {
	namespace score {

		/// \brief TODOCUMENT
		class score_value_reader final {
		private:

		public:
			static aligned_pair_score_value_list read_aligned_pair_score_list_from_property_tree(const boost::property_tree::ptree &);
			static aligned_pair_score_value_list read(std::istream &);
			static aligned_pair_score_value_list read(const boost::filesystem::path &);
		};

	}
}

#endif

