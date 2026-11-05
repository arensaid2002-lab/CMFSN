#pragma once

#include <Kokkos_Core.hpp>
#include <cstddef>

namespace Solver
{
    enum class CellCopyFields
    {
        OldState,
        NewState,
        States,
		velocity,
        Geometry,
        Fluxes,
        All
    };
    //À revoir 
        namespace detail
        {
            // The first dimension is always the cell index.
            template <class ViewType>
            KOKKOS_INLINE_FUNCTION
            void copy_cell_view(const ViewType& view, int k_source, int k_destination)
            {
                if constexpr (ViewType::rank() == 1)
                {
                    view(k_destination) = view(k_source);
                }
                else if constexpr (ViewType::rank() == 2)
                {
                    for (std::size_t face = 0; face < view.extent(1); ++face)
                        view(k_destination, face) = view(k_source, face);
                }
                else if constexpr (ViewType::rank() == 3)
                {
                    for (std::size_t component = 0; component < view.extent(2); ++component)
                        for (std::size_t face = 0; face < view.extent(1); ++face)
                            view(k_destination, face, component) = view(k_source, face, component);
                }
            }

            template <class... Views>
            KOKKOS_INLINE_FUNCTION
            void copy_cell_views(int k_source, int k_destination, const Views&... views)
            {
                (copy_cell_view(views, k_source, k_destination), ...);
            }
        }

        // Call from the execution space that can access the views' memory.
        // Parallel calls require unique destinations and sources that remain unchanged.
        // Face ordering and vector components are copied without transformation.
        template <CellCopyFields Fields = CellCopyFields::States, class CellDataType>
        KOKKOS_INLINE_FUNCTION
        void copy_cell(const CellDataType& cellsdata, int k_source, int k_destination)
        {
            if (k_source == k_destination)
                return;

            if constexpr (Fields == CellCopyFields::OldState || Fields == CellCopyFields::States || Fields == CellCopyFields::All)
            {
                detail::copy_cell_views(k_source, k_destination,
                    cellsdata.rho_old, cellsdata.Momentum_u_old, cellsdata.Momentum_v_old, cellsdata.E_old);
            }
            if constexpr (Fields == CellCopyFields::NewState || Fields == CellCopyFields::States || Fields == CellCopyFields::All)
            {
                detail::copy_cell_views(k_source, k_destination,
                    cellsdata.rho_new, cellsdata.Momentum_u_new, cellsdata.Momentum_v_new, cellsdata.E_new);
            }
			if constexpr (Fields == CellCopyFields::velocity || Fields == CellCopyFields::OldState || Fields == CellCopyFields::NewState || Fields == CellCopyFields::States || Fields == CellCopyFields::All)
			{
				detail::copy_cell_views(k_source, k_destination,
					cellsdata.velocity_u, cellsdata.velocity_v);
			}
            if constexpr (Fields == CellCopyFields::Geometry || Fields == CellCopyFields::All)
            {
                detail::copy_cell_views(k_source, k_destination,
                    cellsdata.Cell_Areas, cellsdata.Normals);
            }
            if constexpr (Fields == CellCopyFields::Fluxes || Fields == CellCopyFields::All)
            {
                detail::copy_cell_views(k_source, k_destination,
                    cellsdata.Flux_mass, cellsdata.Flux_Momentum_u, cellsdata.Flux_Momentum_v, cellsdata.Flux_E, cellsdata.D);
            }
        }
}
