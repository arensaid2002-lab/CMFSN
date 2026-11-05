#pragma once

#include <iostream>
#include "CellsHandler/CellsData.h"
#include "CellOps/CellCopy.h"
#include "NumOps/ArrayOps.h"
#include "NumOps/NumFunctions.h"

using ExecSpace = Kokkos::DefaultExecutionSpace;
using namespace Solver::NumOps;
namespace Solver::CellOps
{
    enum class DataFields
    {
        #pragma region Statefields
        Conservative_Old,
        Conservative_New,
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

        template <typename T, DataFields Fields = DataFields::Velocity>
        KOKKOS_FUNCTION
        auto extract_DataField(CellsData<T> cellsdata, int k, int direction = 0)
        {
			if constexpr (Fields == DataFields::EdgeLength)
			{
				T S = cellsdata.Normals(k, direction, 2);
				return S;
			}
            else if constexpr (Fields == DataFields::Velocity)
            {    
                Kokkos::Array<T, 2> Velocity = { cellsdata.velocity_u(k), cellsdata.velocity_v(k) };
                return Velocity;
            }
            else if constexpr (Fields == DataFields::NormalUnitVector)
            {
                Kokkos::Array<T, 2> Normals = { cellsdata.Normals(k, direction, 0), cellsdata.Normals(k, direction, 1) };
                return Normals;
            }
            else if constexpr (Fields == DataFields::Fluxes)
            {               
                Kokkos::Array<T, 4> Fluxes = { cellsdata.Flux_mass(k, direction), cellsdata.Flux_Momentum_u(k, direction), cellsdata.Flux_Momentum_v(k, direction), cellsdata.Flux_E(k, direction) };
                return Fluxes;
            }
            else if constexpr (Fields == DataFields::Conservative_Old)
            {
                Kokkos::Array<T, 4> Conservative = { cellsdata.rho_old(k), cellsdata.Momentum_u_old(k), cellsdata.Momentum_v_old(k), cellsdata.E_old(k) };
                return Conservative;
            }
            else
            {
                static_assert(Fields == DataFields::EdgeLength || Fields == DataFields::Velocity || Fields == DataFields::NormalUnitVector || Fields == DataFields::Fluxes || Fields == DataFields::Conservative_Old, "Invalid data field specified");
            }
        }
        template <typename T, DataFields Fields = DataFields::Velocity, typename InsertData>
        KOKKOS_FUNCTION
        void Insert_DataField(InsertData data, CellsData<T> cellsdata, int k, int direction = 0)
        {
            if constexpr (Fields == DataFields::EdgeLength && std::is_arithmetic_v<InsertData>)
            {
                cellsdata.Normals(k, direction, 2) = data;
            }
            else if constexpr (Fields == DataFields::Velocity && (std::is_same_v <InsertData, Kokkos::Array<T, 2>>))
            {
				cellsdata.velocity_u(k) = data[0];
				cellsdata.velocity_v(k) = data[1];
            }
            else if constexpr (Fields == DataFields::NormalUnitVector && (std::is_same_v <InsertData, Kokkos::Array<T, 2>>))
            {
                cellsdata.Normals(k, direction, 0) = data[0];
                cellsdata.Normals(k, direction, 1) = data[1];
            }
            else if constexpr (Fields == DataFields::Fluxes && (std::is_same_v <InsertData, Kokkos::Array<T, 4>>))
            {
                cellsdata.Flux_mass(k, direction) = data[0];
                cellsdata.Flux_Momentum_u(k, direction) = data[1];
                cellsdata.Flux_Momentum_v(k, direction) = data[2];
                cellsdata.Flux_E(k, direction) = data[3];
            }
			else if constexpr (Fields == DataFields::Conservative_Old && (std::is_same_v <InsertData, Kokkos::Array<T, 4>>))
			{
				cellsdata.rho_old(k) = data[0];
				cellsdata.Momentum_u_old(k) = data[1];
				cellsdata.Momentum_v_old(k) = data[2];
				cellsdata.E_old(k) = data[3];
			}
			else if constexpr (Fields == DataFields::Conservative_New && (std::is_same_v <InsertData, Kokkos::Array<T, 4>>))
			{
				cellsdata.rho_new(k) = data[0];
				cellsdata.Momentum_u_new(k) = data[1];
				cellsdata.Momentum_v_new(k) = data[2];
				cellsdata.E_new(k) = data[3];
			}
            else
            {
                static_assert((Fields == DataFields::EdgeLength && std::is_arithmetic_v<InsertData>) || (Fields == DataFields::Velocity && (std::is_same_v <InsertData, Kokkos::Array<T, 2>>)) || (Fields == DataFields::NormalUnitVector && (std::is_same_v <InsertData, Kokkos::Array<T, 2>>)) || (Fields == DataFields::Fluxes && (std::is_same_v <InsertData, Kokkos::Array<T, 4>>)) || (Fields == DataFields::Conservative_Old && (std::is_same_v <InsertData, Kokkos::Array<T, 4>>)) || (Fields == DataFields::Conservative_New && (std::is_same_v <InsertData, Kokkos::Array<T, 4>>)), "Invalid data field specified");
            }
        }

        template <typename ExecSpace, CellCopyFields Fields = CellCopyFields::All, typename T>
        void ApplyPeriodicity(CellsData<T> cellsdata)
        {
            Kokkos::parallel_for("Periodicity", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_ghostlayers* cellsdata.Nj_cells), KOKKOS_LAMBDA(int K_g)
            {
                int i_g = K_g % cellsdata.N_ghostlayers;
                int j_g = K_g / cellsdata.N_ghostlayers;

				int i_Ldestination = i_g;
				int j_Ldestination = j_g + cellsdata.N_ghostlayers;
			    int i_Lsource = i_Ldestination + cellsdata.Ni_cells;
                int j_Lsource = j_Ldestination;

                int k_Lcell_source = Cell_index(i_Lsource, j_Lsource, cellsdata.Nitotal);
                int k_Lcell_destination = Cell_index(i_Ldestination, j_Ldestination, cellsdata.Nitotal);
				copy_cell <Fields>(cellsdata, k_Lcell_source, k_Lcell_destination);

                int i_Rsource = i_g + cellsdata.N_ghostlayers;
				int j_Rsource = j_g + cellsdata.N_ghostlayers;
				int i_Rdestination = i_Rsource + cellsdata.Ni_cells;
				int j_Rdestination = j_Rsource;
                
                int k_Rcell_source = Cell_index(i_Rsource, j_Rsource, cellsdata.Nitotal);
                int k_Rcell_destination = Cell_index(i_Rdestination, j_Rdestination, cellsdata.Nitotal);
                copy_cell <Fields>(cellsdata, k_Rcell_source, k_Rcell_destination); 
            });
        }
		
        template<typename ExecSpace, typename T>
        void ApplyInitialConditions(CellsData<T> cellsdata, const initialization_data<T>& init_data)
        {
            Kokkos::parallel_for("InitialConditions", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_cells_total), KOKKOS_LAMBDA(int k)
            {
                cellsdata.rho_old(k) = init_data.rho_inf;
				cellsdata.velocity_u(k) = init_data.velocity_u_inf;
				cellsdata.velocity_v(k) = init_data.velocity_v_inf;
                cellsdata.Momentum_u_old(k) = init_data.rho_inf * init_data.velocity_u_inf;
                cellsdata.Momentum_v_old(k) = init_data.rho_inf * init_data.velocity_v_inf;
                cellsdata.E_old(k) = (T)1 / (init_data.gamma_ - (T)1) + (T)0.5 * init_data.gamma_*init_data.Mach_inf_*init_data.Mach_inf_;
            });
        }

        template<typename ExecSpace, typename T>
        void ApplyNonPenetrationConditions(CellsData<T> cellsdata, const initialization_data<T>& init_data)
        {
            Kokkos::parallel_for("NonPenetrationConditions", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.Ni_cells * cellsdata.N_ghostlayers), KOKKOS_LAMBDA(int k_g)
            {
                int i_g = k_g % cellsdata.Ni_cells;
                int j_g = k_g / cellsdata.Ni_cells;

                int i_source = i_g + cellsdata.N_ghostlayers;
                int j_source = 2 * cellsdata.N_ghostlayers - 1 - j_g;

                int i_destination = i_source;
                int& j_destination = j_g;


                int k_source = Cell_index(i_source, j_source, cellsdata.Nitotal);
                int k_destination = Cell_index(i_destination, j_destination, cellsdata.Nitotal);
				int k_wall = Cell_index(i_destination, cellsdata.N_ghostlayers, cellsdata.Nitotal);
                int& k_ghost = k_destination;
                
				//copy the old state from the source cell to the ghost cell
                copy_cell <CellCopyFields::OldState>(cellsdata, k_source, k_destination);

                //eval velocity reflection
                auto ghost_normal = -1 * extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k_wall, 0);
                auto RealCell_velocity = extract_DataField<T, DataFields::Velocity>(cellsdata, k_source);
                auto GhostCell_velocity = RealCell_velocity - (T)2 * dot(RealCell_velocity, ghost_normal) * ghost_normal;
                Insert_DataField<T, DataFields::Velocity>(GhostCell_velocity, cellsdata, k_destination);
				
				//update state and fluxes in the ghost cell based on the reflected velocity
                T pressure = eval_pressure<T>(cellsdata.rho_old(k_ghost), cellsdata.E_old(k_ghost), cellsdata.velocity_u(k_ghost), cellsdata.velocity_v(k_ghost), init_data.gamma_);
				cellsdata.Momentum_u_old(k_ghost) = cellsdata.rho_old(k_ghost) * cellsdata.velocity_u(k_ghost);
				cellsdata.Momentum_v_old(k_ghost) = cellsdata.rho_old(k_ghost) * cellsdata.velocity_v(k_ghost);
				cellsdata.Flux_Momentum_u(k_ghost, 2) = ghost_normal[0] * pressure;
				cellsdata.Flux_Momentum_v(k_ghost, 2) = ghost_normal[1] * pressure;
            });
        }
        
        template <typename ExecSpace, typename T>
        void ApplyRiemanInvariantConditions(CellsData<T> cellsdata, const initialization_data<T>& init_data)
        {
            if (cellsdata.N_ghostlayers <= 0)
                return;

            Kokkos::parallel_for("RiemanInvariant", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.Ni_cells), KOKKOS_LAMBDA(int i_g)
            {
                int i_cells = i_g + cellsdata.N_ghostlayers;
                int j_Cell = cellsdata.N_ghostlayers + cellsdata.Nj_cells - 1;
                int k_cell = Cell_index(i_cells, j_Cell, cellsdata.Nitotal);
                
                int j_ghost_first = j_Cell + 1;
                int j_ghost_last = j_ghost_first + 1;
                int k_ghost_first = Cell_index(i_cells, j_ghost_first, cellsdata.Nitotal);
                int k_ghost_last = Cell_index(i_cells, j_ghost_last, cellsdata.Nitotal);

                auto normal = extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k_cell, 2);
                auto V_in = extract_DataField<T, DataFields::Velocity>(cellsdata, k_cell);
                auto V_inf = Kokkos::Array<T, 2>{ init_data.velocity_u_inf, init_data.velocity_v_inf };

                T P_in = eval_pressure<T>(cellsdata.rho_old(k_cell), cellsdata.E_old(k_cell), V_in[0], V_in[1], init_data.gamma_);
                T c_in = eval_sound_speed<T>(init_data.gamma_, P_in, cellsdata.rho_old(k_cell));
                auto boundary_state = eval_RiemannBoundaryState<T>(V_in, V_inf, normal, c_in, init_data.c_inf, init_data.gamma_, cellsdata.rho_old(k_cell), init_data.rho_inf);

                auto V_boundary = boundary_state.first;
                T rho_boundary = boundary_state.second.first;
                T pressure_boundary = boundary_state.second.second;
                T Energy = eval_Energy<T>(rho_boundary, V_boundary[0], V_boundary[1], pressure_boundary, init_data.gamma_);

                auto Cell_State = extract_DataField<T, DataFields::Conservative_Old>(cellsdata, k_cell);
                Kokkos::Array<T, 4> Boundary_States = {rho_boundary, rho_boundary * V_boundary[0], rho_boundary * V_boundary[1], Energy};
                auto Ghost_State_First = T{2} * Boundary_States - Cell_State;
                auto Ghost_State_Last = Ghost_State_First;
                Kokkos::Array<T, 2> Ghost_Velocity = { Ghost_State_First[1] / Ghost_State_First[0], Ghost_State_First[2] / Ghost_State_First[0] };

                Insert_DataField<T, DataFields::Conservative_Old>(Ghost_State_First, cellsdata, k_ghost_first);
                Insert_DataField<T, DataFields::Velocity>(Ghost_Velocity, cellsdata, k_ghost_first);
                if (cellsdata.N_ghostlayers > 1)
                {
                    Insert_DataField<T, DataFields::Conservative_Old>(Ghost_State_Last, cellsdata, k_ghost_last);
                    Insert_DataField<T, DataFields::Velocity>(Ghost_Velocity, cellsdata, k_ghost_last);
                }
                
            });
        }



        /**
        * Evaluates the geometrical data for a given cell.
        * @param i The cell index in the i-direction.
        * @param j The cell index in the j-direction.
        * @param direction The direction of the normal vector.
        * @param meshdata The mesh data.
        */
        template<typename T>
        KOKKOS_FUNCTION
        Kokkos::Array<T, 3> eval_CellNormal(int i, int j, int direction, const MeshData_Structured<T>& meshdata)
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
			Kokkos::Array<T, 3> Normal = { delta_y / Norm, -delta_x / Norm, Norm };
            return Normal;
        }



        template <typename ExecSpace, typename T>
        void eval_Residual()
        {

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
                int i = k_origin % cellsdata.Ni_cells;
                int j = k_origin / cellsdata.Ni_cells;

                int i_cell = i + cellsdata.N_ghostlayers;
                int j_cell = j + cellsdata.N_ghostlayers;
                int k_cell = Cell_index(i_cell, j_cell, cellsdata.Nitotal);
                cellsdata.Cell_Areas(k_cell) = eval_CellArea(i, j, meshdata);

                for (int direction = 0; direction < cellsdata.N_edges; direction++)
                {
                    const auto normals = eval_CellNormal(i, j, direction, meshdata);

                    cellsdata.Normals(k_cell, direction, 0) = normals[0];
                    cellsdata.Normals(k_cell, direction, 1) = normals[1];
                    cellsdata.Normals(k_cell, direction, 2) = normals[2];
                }
            });
        }
        


}
