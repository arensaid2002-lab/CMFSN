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
                Kokkos::Array<T, 2> Velocity = {
                    cellsdata.Momentum_u_old(k) / cellsdata.rho_old(k),
                    cellsdata.Momentum_v_old(k) / cellsdata.rho_old(k)
                };
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
        void Update_DataField(InsertData data, CellsData<T> cellsdata, int k, int direction = 0)
        {
            if constexpr (Fields == DataFields::EdgeLength && std::is_arithmetic_v<InsertData>)
            {
                cellsdata.Normals(k, direction, 2) = data;
            }
            else if constexpr (Fields == DataFields::Velocity && (std::is_same_v <InsertData, Kokkos::Array<T, 2>>))
            {
                //Here Energy update is omitted, it is assumed the velocity magnitude remains constant
                cellsdata.velocity_u(k) = data[0];
                cellsdata.velocity_v(k) = data[1];

                cellsdata.Momentum_u_old(k) = cellsdata.rho_old(k) * data[0];
                cellsdata.Momentum_v_old(k) = cellsdata.rho_old(k) * data[1]; 
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

                const int k_Lcell_destination = Cell_index(i_g, j_g + cellsdata.N_ghostlayers, cellsdata.Nitotal);
                const int k_Lcell_source = k_Lcell_destination + cellsdata.Ni_cells;
                const int k_Rcell_source = k_Lcell_destination + cellsdata.N_ghostlayers;
                const int k_Rcell_destination = k_Rcell_source + cellsdata.Ni_cells;

                copy_cell <Fields>(cellsdata, k_Lcell_source, k_Lcell_destination);
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

                const int stride = cellsdata.Nitotal - 1;
                const int k_wall = Cell_index(i_g + cellsdata.N_ghostlayers, cellsdata.N_ghostlayers, cellsdata.Nitotal);
                const int k_source = k_wall + (cellsdata.N_ghostlayers - 1 - j_g) * stride;
                const int k_destination = k_wall - (cellsdata.N_ghostlayers - j_g) * stride;
                const int k_ghost = k_destination;
                
				//copy the old state from the source cell to the ghost cell
                copy_cell <CellCopyFields::OldState>(cellsdata, k_source, k_destination);

                //eval velocity reflection
                auto ghost_normal = -1 * extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k_wall, 0);
                auto RealCell_velocity = extract_DataField<T, DataFields::Velocity>(cellsdata, k_source);
                auto GhostCell_velocity = RealCell_velocity - (T)2 * dot(RealCell_velocity, ghost_normal) * ghost_normal;
                Update_DataField<T, DataFields::Velocity>(GhostCell_velocity, cellsdata, k_ghost);
                update_Stateflux<T>(cellsdata, init_data)
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

                const int stride = cellsdata.Nitotal - 1;
                const int k_cell = Cell_index(i_cells, j_Cell, cellsdata.Nitotal);
                const int k_ghost_first = k_cell + stride;
                const int k_ghost_last = k_cell + 2 * stride;

                auto normal = extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k_cell, 2);
                auto V_in = extract_DataField<T, DataFields::Velocity>(cellsdata, k_cell);
                auto V_inf = Kokkos::Array<T, 2>{ init_data.velocity_u_inf, init_data.velocity_v_inf };

                T P_in = cellsdata.Pressure(k_cell);
                T c_in = cellsdata.SoundSpeed(k_cell);
                auto boundary_state = eval_RiemannBoundaryState<T>(V_in, V_inf, normal, c_in, init_data.c_inf, init_data.gamma_, cellsdata.rho_old(k_cell), init_data.rho_inf);

                auto V_boundary = boundary_state.first;
                auto rho_pressure = boundary_state.second;
                T rho_boundary = rho_pressure.first;
                T pressure_boundary = rho_pressure.second;
                T Energy = eval_Energy<T>(rho_boundary, V_boundary[0], V_boundary[1], pressure_boundary, init_data.gamma_);

                auto Cell_State = extract_DataField<T, DataFields::Conservative_Old>(cellsdata, k_cell);
                Kokkos::Array<T, 4> Boundary_States = {rho_boundary, rho_boundary * V_boundary[0], rho_boundary * V_boundary[1], Energy};
                auto Ghost_State_First = T{2} * Boundary_States - Cell_State;
                auto Ghost_State_Last = Ghost_State_First;
                Kokkos::Array<T, 2> Ghost_Velocity = { Ghost_State_First[1] / Ghost_State_First[0], Ghost_State_First[2] / Ghost_State_First[0] };

                Update_DataField<T, DataFields::Conservative_Old>(Ghost_State_First, cellsdata, k_ghost_first);
                Update_DataField<T, DataFields::Conservative_Old>(Ghost_State_Last, cellsdata, k_ghost_last);
                
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
            Kokkos::parallel_for("eval_geometrical_data", Kokkos::RangePolicy<ExecSpace>(0, meshdata.N_Cells), KOKKOS_LAMBDA(int k_mesh)
            {
                int i = k_mesh % cellsdata.Ni_cells;
                int j = k_mesh / cellsdata.Ni_cells;

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
        
        template <typename ExecSpace, typename T>
        void eval_DerivedFields(CellsData<T> cellsdata, const initialization_data<T>& init_data)
        {
            Kokkos::parallel_for("eval_DerivedFields", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_cells_mesh), KOKKOS_LAMBDA(int k_mesh)
            {
                int i = k_mesh % cellsdata.Ni_cells;
                int j = k_mesh / cellsdata.Ni_cells;

                int i_cell = i + cellsdata.N_ghostlayers;
                int j_cell = j + cellsdata.N_ghostlayers;
                int k_cell = Cell_index(i_cell, j_cell, cellsdata.Nitotal);
                
                T Pressure_value = eval_pressure<T>(cellsdata.rho_old(k_cell), cellsdata.E_old(k_cell), cellsdata.velocity_u(k_cell), cellsdata.velocity_v(k_cell), init_data.gamma_);
                T Sound_speed_value = eval_sound_speed<T>(init_data.gamma_, Pressure_value, cellsdata.rho_old(k_cell));

                cellsdata.Pressure(k_cell) = Pressure_value;
                cellsdata.SoundSpeed(k_cell) = Sound_speed_value;
                for (int direction = 0; direction < cellsdata.N_edges; direction++)
                {
                    auto normal = extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k_cell, direction);
                    auto Velocity = extract_DataField<T, DataFields::Velocity>(cellsdata, k_cell, direction);
                    T EdgeLength = extract_DataField<T, DataFields::EdgeLength>(cellsdata, k_cell, direction);
                    T Lambda = abs(dot(Velocity, normal)) * EdgeLength;
                    cellsdata.Lambdas(k_cell, direction) = Lambda;
                }
            });
        }


        template <typename T>
        KOKKOS_INLINE_FUNCTION
        void update_Stateflux(CellsData<T> cellsdata, const initialization_data<T>& init_data,int k_cell, int direction)
        {
            auto velocity = extract_DataField<T, DataFields::Velocity>(cellsdata, k_cell);
            auto normal = extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k_cell, direction);
            auto Conservative_state = extract_DataField<T, DataFields::Conservative_Old>(cellsdata, k_cell);
            T Pressure = cellsdata.Pressure(k_cell);
            T Vn = dot(velocity, normal);

            Kokkos::Array<T, 4> Pressure_Energy = { 0, normal[0] * Pressure, normal[1] * Pressure, Pressure * Vn };
            auto flux = Conservative_state * Vn + Pressure_Energy;
            Update_DataField<T, DataFields::Fluxes>(flux, cellsdata, k_cell, direction);
        }
        
        template <typename T>
        KOKKOS_INLINE_FUNCTION
        Kokkos::Array<T, 4> eval_Dissipation(CellsData<T> cellsdata, const initialization_data<T>& init_data, int k_cell)
        {
            const int stride = cellsdata.Nitotal - 1;
            const int k_cell_0 = k_cell - stride;
            const int k_cell_1 = k_cell + 1;
            const int k_cell_2 = k_cell + stride;
            const int k_cell_3 = k_cell - 1;

           
            Kokkos::Array<T, 4> Lambdas_2; 
            Lambdas_2[0] = mean<T>(cellsdata.Lambdas(k_cell, 0), cellsdata.Lambdas(k_cell_0, 2));
            Lambdas_2[1] = mean<T>(cellsdata.Lambdas(k_cell, 1), cellsdata.Lambdas(k_cell_1, 3));
            Lambdas_2[2] = mean<T>(cellsdata.Lambdas(k_cell, 2), cellsdata.Lambdas(k_cell_2, 0));
            Lambdas_2[3] = mean<T>(cellsdata.Lambdas(k_cell, 3), cellsdata.Lambdas(k_cell_3, 1));



            T pressure_captor_i = cellsdata.Pressure(k_cell_0) + 
            
        }



}
