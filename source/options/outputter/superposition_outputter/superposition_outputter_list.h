/// \file
/// \brief The superposition_outputter_list class header

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

#ifndef SUPERPOSITION_OUTPUTTER_LIST_H_INCLUDED
#define SUPERPOSITION_OUTPUTTER_LIST_H_INCLUDED

#include <boost/ptr_container/ptr_vector.hpp>

#include <iosfwd>

namespace cath { namespace sup { class superposition_context; } }
namespace cath { namespace opts { class superposition_outputter; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class superposition_outputter_list final {
		private:
			boost::ptr_vector<superposition_outputter> outputters;

		public:
			void push_back(const superposition_outputter &);
			bool empty() const;

			// Provide iterators
			using const_iterator = boost::ptr_vector<superposition_outputter>::const_iterator;
			const_iterator begin() const;
			const_iterator end() const;
		};

		void use_all_superposition_outputters(const superposition_outputter_list &,
		                                      const sup::superposition_context &,
		                                      std::ostream &,
		                                      std::ostream &);

		bool any_superposition_outputters_involve_display_spec(const superposition_outputter_list &);
	}
}

#endif
