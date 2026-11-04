#pragma once

#include <iostream>
#include "CellsHandler/CellsData.h"
#include "CellOps/CellCopy.h"

using ExecSpace = Kokkos::DefaultExecutionSpace;
namespace Solver::CellOps
{
    enum class DataFields
    {
        #pragma region Statefields
        Conservative,
        Velocity,
        Fluxes,
        #pragma endregion
        #pragma region Geometrical Data Fields
        NormalUnitVector,
        EdgeLength
        #pragma endregion
    };

        KOKKOS_INLINE_FUNCTION
        int Node_index(int i, int j, int Ni) {return i + Ni*j;}

        KOKKOS_INLINE_FUNCTION
        int Cell_index(int i, int j, int Ni) {return i + (Ni-1)*j;}

        template <typename T, std::size_t N, DataFields Fields = DataFields::Velocity>
        KOKKOS_FUNCTION
        Kokkos::Array<T, N> extract_DataField_toArray(CellsData<T> cellsdata, int k, int direction = 0)
        {
            if constexpr (Fields == DataFields::Velocity)
            {
                N = 2; // size of velocity vector
                Kokkos::Array<T, N> Velocity = {cellsdata.velocity_u(k), cellsdata.velocity_v(k)};
                return Velocity;
            }
            else if constexpr (Fields == DataFields::NormalUnitVector)
            {
                N = 2; // size of normal vector
                Kokkos::Array<T, N> Normals = {cellsdata.Normals(k, direction, 0), cellsdata.Normals(k, direction, 1)};
                return Normals;
            }
            else if constexpr (Fields == DataFields::Fluxes)
            {
                N = 4; // size of flux vector
                Kokkos::Array<T, N> Fluxes = {cellsdata.Flux_mass(k, direction), cellsdata.Flux_Momentum_u(k, direction), cellsdata.Flux_Momentum_v(k, direction), cellsdata.Flux_E(k, direction)};
                return Fluxes;
            }
            else if constexpr (Fields == DataFields::Conservative)
            {
                N = 4; // size of Conservative vector
                Kokkos::Array<T, N> Conservative = {cellsdata.rho_old(k), cellsdata.Momentum_u_old(k), cellsdata.Momentum_v_old(k), cellsdata.E_old(k) };
                return Conservative;
            }
            else
            {
                static_assert(Fields == DataFields::Velocity|| Fields == DataFields::NormalUnitVector || Fields == DataFields::Fluxes || Fields == DataFields::Conservative, "Invalid data field specified");
            }
        }

        template <typename ExecSpace, CellCopyFields Fields = CellCopyFields::All, typename T>
        void ApplyPeriodicity(CellsData<T> cellsdata, const MeshData_Structured<T>& meshdata)
        {
            Kokkos::parallel_for("Periodicity", Kokkos::RangePolicy<ExecSpace>(0, meshdata.Nj-1), KOKKOS_LAMBDA(int j)
            {

			    int j_cell = j + cellsdata.N_ghostlayers;
                for (int i_cell_destination = 0; i_cell_destination < cellsdata.N_ghostlayers; i_cell_destination++)
                {
                    int i_source = (meshdata.Ni - 1) + i_cell_destination;
				    int k_cell_source = Cell_index(i_source, j_cell, cellsdata.Nitotal);
                    int k_cell_destination = Cell_index(i_cell_destination, j_cell, cellsdata.Nitotal);
                    copy_cell <Fields>(cellsdata, k_cell_source, k_cell_destination);
                }
                for (int i_cell_destination = (cellsdata.Nitotal-1)-cellsdata.N_ghostlayers; i_cell_destination < cellsdata.Nitotal-1; i_cell_destination++)
                {
                    int i_source = i_cell_destination - (meshdata.Ni - 1);
                    int k_cell_source = Cell_index(i_source, j_cell, cellsdata.Nitotal);
                    int k_cell_destination = Cell_index(i_cell_destination, j_cell, cellsdata.Nitotal);
                    copy_cell <Fields>(cellsdata, k_cell_source, k_cell_destination);
                }
            });
            Kokkos::fence();
        }
		
        template<typename ExecSpace, typename T>
        void ApplyInitialConditions(CellsData<T> cellsdata, const initialization_data<T>& init_data)
        {
            Kokkos::parallel_for("InitialConditions", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_cells_total), KOKKOS_LAMBDA(int k)
            {
                cellsdata.rho_old(k) = init_data.rho_inf_;
                cellsdata.Momentum_u_old(k) = init_data.rho_inf_ * init_data.velocity_u_inf_;
                cellsdata.Momentum_v_old(k) = init_data.rho_inf_ * init_data.velocity_v_inf_;
                cellsdata.E_old(k) = 1.0f / (init_data.gamma_ - 1.0f) + 0.5f * init_data.gamma_*init_data.Mach_inf_*init_data.Mach_inf_;
            });
            Kokkos::fence();
        }

        template<typename ExecSpace, typename T>
        void ApplyPenetrationConditions(CellsData<T> cellsdata);

        /**
        * Evaluates the geometrical data for a given cell.
        * @param i The cell index in the i-direction.
        * @param j The cell index in the j-direction.
        * @param direction The direction of the normal vector.
        * @param meshdata The mesh data.
        */
        template<typename T, DataFields Field = DataFields::NormalUnitVector>
        KOKKOS_FUNCTION
        auto eval_CellNormal(int i, int j, int direction, const MeshData_Structured<T>& meshdata)
        {
            int Ni = meshdata.Ni;
            T delta_x = T{ 0 }, delta_y = T{ 0 };

            switch (direction)
            {
            case 0: // Bottom
                delta_x = meshdata.X(Node_index(i + 1, j, Ni)) - meshdata.X(Node_index(i, j, Ni));
                delta_y = meshdata.Y(Node_index(i + 1, j, Ni)) - meshdata.Y(Node_index(i, j, Ni));

                break;

            case 1: // Right
                delta_x = meshdata.X(Node_index(i + 1, j + 1, Ni)) - meshdata.X(Node_index(i + 1, j, Ni));
                delta_y = meshdata.Y(Node_index(i + 1, j + 1, Ni)) - meshdata.Y(Node_index(i + 1, j, Ni));

                break;
            case 2: // Top
                delta_x = meshdata.X(Node_index(i, j + 1, Ni)) - meshdata.X(Node_index(i + 1, j + 1, Ni));
                delta_y = meshdata.Y(Node_index(i, j + 1, Ni)) - meshdata.Y(Node_index(i + 1, j + 1, Ni));

                break;
            case 3: // Left
                delta_x = meshdata.X(Node_index(i, j, Ni)) - meshdata.X(Node_index(i, j + 1, Ni));
                delta_y = meshdata.Y(Node_index(i, j, Ni)) - meshdata.Y(Node_index(i, j + 1, Ni));

                break;
            }
            T Norm = Kokkos::sqrt(delta_x * delta_x + delta_y * delta_y);
            if constexpr (Field == DataFields::NormalUnitVector)
            {
                return Kokkos::pair<T,T> { delta_y / Norm, -delta_x / Norm };
            }
            else if constexpr (Field == DataFields::EdgeLength)
            {
                return Norm;
            }
            else
            {
                static_assert(Field == DataFields::NormalUnitVector || Field == DataFields::EdgeLength, "Invalid GeometricalDataFields value");
            }
        }

        /**
        * Evaluates the area of a given cell.
        * @param i The cell index in the i-direction.
        * @param j The cell index in the j-direction.
        * @param meshdata The mesh data.
        * @return The area of the cell.
        */
        template <typename T>
        KOKKOS_FUNCTION
        T eval_CellArea(int i, int j, const MeshData_Structured<T>& meshdata)
        {
            int Ni = meshdata.Ni;
            T delta_x1 = meshdata.X(Node_index(i + 1, j + 1, Ni)) - meshdata.X(Node_index(i, j, Ni));
            T delta_y1 = meshdata.Y(Node_index(i + 1, j + 1, Ni)) - meshdata.Y(Node_index(i, j, Ni));

            T delta_x2 = meshdata.X(Node_index(i, j + 1, Ni)) - meshdata.X(Node_index(i + 1, j, Ni));
            T delta_y2 = meshdata.Y(Node_index(i, j + 1, Ni)) - meshdata.Y(Node_index(i + 1, j, Ni));

            return T{ 0.5 } * ((delta_x1 * delta_y2 - delta_x2 * delta_y1));
        }

        template <typename ExecSpace, typename T>
        void eval_GeometricalData(CellsData<T> cellsdata, const MeshData_Structured<T>& meshdata)
        {
            Kokkos::parallel_for("eval_geometrical_data", Kokkos::RangePolicy<ExecSpace>(0, meshdata.N_Cells), KOKKOS_LAMBDA(int k_origin)
            {
                int i = k_origin % (meshdata.Ni - 1);
                int j = k_origin / (meshdata.Ni - 1);

                int i_cell = i + cellsdata.N_gCell_layers;
                int j_cell = j + cellsdata.N_gCell_layers;
                int k_cell = Cell_index(i_cell, j_cell, cellsdata.Nitotal);
                cellsdata.Cell_Areas(k_cell) = eval_Cell_Area(i, j, meshdata);


                for (int direction = 0; direction < meshdata.N_edges; direction++)
                {
                    const auto normals = eval_Cell_Normal(i, j, direction, meshdata);

                    cellsdata.Normals(k_cell, direction, 0) = normals[0];
                    cellsdata.Normals(k_cell, direction, 1) = normals[1];
                }
            });
            Kokkos::fence();
        }
        
        template <typename T, std::size_t N, DataFields Fields = DataFields::Velocity>
        KOKKOS_FUNCTION
        Kokkos::Array<T, N> extract_DataField_toArray(CellsData<T> cellsdata, int k, int direction = 0)
        {
            if constexpr (Fields == DataFields::Velocity)
            {
                N = 2; // size of velocity vector
                Kokkos::Array<T, N> Velocity = { cellsdata.velocity_u(k), cellsdata.velocity_v(k) };
                return Velocity;
            }
            else if constexpr (Fields == DataFields::NormalUnitVector)
            {
                N = 2; // size of normal vector
                Kokkos::Array<T, N> Normals = { cellsdata.Normals(k, direction, 0), cellsdata.Normals(k, direction, 1) };
                return Normals;
            }
            else if constexpr (Fields == DataFields::Fluxes)
            {
                N = 4; // size of flux vector
                Kokkos::Array<T, N> Fluxes = { cellsdata.Flux_mass(k, direction), cellsdata.Flux_Momentum_u(k, direction), cellsdata.Flux_Momentum_v(k, direction), cellsdata.Flux_E(k, direction) };
                return Fluxes;
            }
            else if constexpr (Fields == DataFields::Conservative)
            {
                N = 4; // size of Conservative vector
                Kokkos::Array<T, N> Conservative = { cellsdata.rho_old(k), cellsdata.Momentum_u_old(k), cellsdata.Momentum_v_old(k), cellsdata.E_old(k) };
                return Conservative;
            }
            else
            {
                static_assert(Fields == DataFields::Velocity || Fields == DataFields::NormalUnitVector || Fields == DataFields::Fluxes || Fields == DataFields::Conservative, "Invalid data field specified");
            }
        }

}
