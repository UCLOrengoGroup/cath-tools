/// \file
/// \brief The maybe_unused_namespace_scope_constexpr header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_DETAIL_MAYBE_UNUSED_NAMESPACE_SCOPE_CONSTEXPR_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_DETAIL_MAYBE_UNUSED_NAMESPACE_SCOPE_CONSTEXPR_HPP

#define MAYBE_UNUSED_NAMESPACE_SCOPE_CONSTEXPR( name )                         \
   namespace pseudo_use__ns__##name {                                          \
      template <typename> struct pseudo_use__struct__##name                {}; \
      template <        > struct pseudo_use__struct__##name<decltype(name)>{}; \
   }

#endif
