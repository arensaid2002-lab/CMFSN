#pragma once

#include <iostream>
#include <stdexcept>
#include <utility>
#include <limits>
#include <cmath>
#include "CellsHandler/CellsData.h"
#include "CellOps/CellCopy.h"
#include "SolverTypes.h"
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
        Temperature,
		DerivedFields,
        Fluxes,
        Dissipation_History,

        #pragma endregion
        #pragma region Geometrical Data Fields
        NormalUnitVector,
        EdgeLength,
        CellCenter
        #pragma endregion
    };

        KOKKOS_INLINE_FUNCTION
        int Node_index(int i, int j, int Ni) {return i + Ni*j;}

        KOKKOS_INLINE_FUNCTION
        int Cell_index(int i, int j, int Ni_total_cells) {return i + Ni_total_cells*j;}

        KOKKOS_INLINE_FUNCTION
        int Cell_index_sides(int k_cell, int row_stride, int direction)
        {
            switch (direction)
            {
            case 0: return k_cell - row_stride;
            case 1: return k_cell + 1;
            case 2: return k_cell + row_stride;
            case 3: return k_cell - 1;
            default: return -1; // Invalid direction
            }
        };



        template <typename T, DataFields Fields = DataFields::Velocity>
        KOKKOS_FUNCTION
        auto extract_DataField(const CellsData<T>& cellsdata, int k, int direction = 0)
        {
			if constexpr (Fields == DataFields::EdgeLength)
			{
				T S = cellsdata.Normals(k, direction, 2);
				return S;
			}
            else if constexpr (Fields == DataFields::Velocity)
            {    
                Kokkos::Array<T, 2> Velocity = {cellsdata.Momentum_u_old(k) / cellsdata.rho_old(k), cellsdata.Momentum_v_old(k) / cellsdata.rho_old(k)
                };
                return Velocity;
            }
            else if constexpr (Fields == DataFields::Temperature)
            {
                return cellsdata.Pressure(k) / cellsdata.rho_old(k);
            }
            else if constexpr (Fields == DataFields::NormalUnitVector)
            {
                Kokkos::Array<T, 2> Normals = { cellsdata.Normals(k, direction, 0), cellsdata.Normals(k, direction, 1) };
                return Normals;
            }
            else if constexpr (Fields == DataFields::CellCenter)
            {
                Kokkos::Array<T, 2> CellCenter = { cellsdata.Cell_Centers(k, 0), cellsdata.Cell_Centers(k, 1) };
                return CellCenter;
            }
            else if constexpr (Fields == DataFields::Fluxes)
            {               
                Kokkos::Array<T, 4> Fluxes = { cellsdata.Flux_mass(k, direction), cellsdata.Flux_Momentum_u(k, direction), cellsdata.Flux_Momentum_v(k, direction), cellsdata.Flux_E(k, direction) };
                return Fluxes;
            }
            else if constexpr (Fields == DataFields::Dissipation_History)
            {
                Kokkos::Array<T, 4> Dissipation = { cellsdata.D_history(k, 0), cellsdata.D_history(k, 1), cellsdata.D_history(k, 2), cellsdata.D_history(k, 3) };
                return Dissipation;
            }
            else if constexpr (Fields == DataFields::Conservative_Old)
            {
                Kokkos::Array<T, 4> Conservative = { cellsdata.rho_old(k), cellsdata.Momentum_u_old(k), cellsdata.Momentum_v_old(k), cellsdata.E_old(k) };
                return Conservative;
            }
            else
            {
                static_assert(Fields == DataFields::EdgeLength || Fields == DataFields::Velocity || Fields == DataFields::Temperature || Fields == DataFields::NormalUnitVector || Fields == DataFields::CellCenter || Fields == DataFields::Fluxes || Fields == DataFields::Dissipation_History || Fields == DataFields::Conservative_Old, "Invalid data field specified");
            }
        }
        
        template <typename T, DataFields Fields>
        KOKKOS_FUNCTION
        void Update_DataField(Kokkos::Array<T, 4> data, const CellsData<T>& cellsdata, int k)
        {
            static_assert(Fields == DataFields::Dissipation_History, "This overload inserts the retained dissipation for one cell");
            cellsdata.D_history(k, 0) = data[0];
            cellsdata.D_history(k, 1) = data[1];
            cellsdata.D_history(k, 2) = data[2];
            cellsdata.D_history(k, 3) = data[3];
        }

        template <typename T, DataFields Fields, typename InsertData>
        KOKKOS_FUNCTION
        void Update_DataField(InsertData data, const CellsData<T>& cellsdata, int k, int direction)
        {
            if constexpr (Fields == DataFields::EdgeLength && std::is_arithmetic_v<InsertData>)
            {
                cellsdata.Normals(k, direction, 2) = data;
            }
            else if constexpr (Fields == DataFields::NormalUnitVector && std::is_same_v<InsertData, Kokkos::Array<T, 2>>)
            {
                cellsdata.Normals(k, direction, 0) = data[0];
                cellsdata.Normals(k, direction, 1) = data[1];
            }
            else if constexpr (Fields == DataFields::Fluxes && std::is_same_v<InsertData, Kokkos::Array<T, 4>>)
            {
                cellsdata.Flux_mass(k, direction) = data[0];
                cellsdata.Flux_Momentum_u(k, direction) = data[1];
                cellsdata.Flux_Momentum_v(k, direction) = data[2];
                cellsdata.Flux_E(k, direction) = data[3];
            }
            else
            {
                static_assert((Fields == DataFields::EdgeLength && std::is_arithmetic_v<InsertData>) || (Fields == DataFields::NormalUnitVector && std::is_same_v<InsertData, Kokkos::Array<T, 2>>) || (Fields == DataFields::Fluxes && std::is_same_v<InsertData, Kokkos::Array<T, 4>>), "This overload accepts face geometry or precomputed fluxes only");
            }
        }
        
        // Refresh the current state's derived fields without inserting a new state.
        template <typename T, DataFields Fields>
        KOKKOS_INLINE_FUNCTION
        void Update_DataField(const CellsData<T>& cellsdata, const initialization_data<T>& init_data, int k, int direction)
        {
            static_assert(Fields == DataFields::DerivedFields, "This overload refreshes derived fields for one face");
            const auto velocity = extract_DataField<T, DataFields::Velocity>(cellsdata, k);
            const T pressure = eval_pressure<T>(cellsdata.rho_old(k), cellsdata.E_old(k), velocity[0], velocity[1], init_data.gamma);
            const T sound_speed = eval_sound_speed<T>(init_data.gamma, pressure, cellsdata.rho_old(k));
            cellsdata.velocity_u(k) = velocity[0];
            cellsdata.velocity_v(k) = velocity[1];
            cellsdata.Pressure(k) = pressure;
            cellsdata.SoundSpeed(k) = sound_speed;
            const auto normal = extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k, direction);
            const T edge_length = extract_DataField<T, DataFields::EdgeLength>(cellsdata, k, direction);
            cellsdata.Lambdas(k, direction) = eval_Lambda<T>(velocity, normal, sound_speed, edge_length);
        }

        template <typename T, DataFields Fields>
        KOKKOS_INLINE_FUNCTION
        void Update_DataField(const CellsData<T>& cellsdata, const initialization_data<T>& init_data, int k)
        {
            static_assert(Fields == DataFields::DerivedFields, "This overload refreshes derived fields for all faces");
            Update_DataField<T, DataFields::DerivedFields>(cellsdata, init_data, k, 0);
            const auto velocity = extract_DataField<T, DataFields::Velocity>(cellsdata, k);
            for (int direction = 1; direction < cellsdata.N_faces; ++direction)
            {
                const auto normal = extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k, direction);
                const T edge_length = extract_DataField<T, DataFields::EdgeLength>(cellsdata, k, direction);
                cellsdata.Lambdas(k, direction) = eval_Lambda<T>(velocity, normal, cellsdata.SoundSpeed(k), edge_length);
            }
        }

        template <typename T, DataFields Fields, typename InsertData>
        KOKKOS_FUNCTION
        void Update_DataField(InsertData data, const CellsData<T>& cellsdata, const initialization_data<T>& init_data, int k)
        {
            if constexpr (Fields == DataFields::Conservative_Old && std::is_same_v<InsertData, Kokkos::Array<T, 4>>)
            {
                cellsdata.rho_old(k) = data[0];
                cellsdata.Momentum_u_old(k) = data[1];
                cellsdata.Momentum_v_old(k) = data[2];
                cellsdata.E_old(k) = data[3];
                Update_DataField<T, DataFields::DerivedFields>(cellsdata, init_data, k);
            }
            else if constexpr (Fields == DataFields::Conservative_New && std::is_same_v<InsertData, Kokkos::Array<T, 4>>)
            {
                cellsdata.rho_new(k) = data[0];
                cellsdata.Momentum_u_new(k) = data[1];
                cellsdata.Momentum_v_new(k) = data[2];
                cellsdata.E_new(k) = data[3];
            }
            else
            {
                static_assert((Fields == DataFields::Conservative_Old || Fields == DataFields::Conservative_New) && std::is_same_v<InsertData, Kokkos::Array<T, 4>>, "This field requires a direction or a different data type");
            }
        }

        template <typename T, DataFields Fields, typename InsertData>
        KOKKOS_FUNCTION
        void Update_DataField(InsertData data, const CellsData<T>& cellsdata, const initialization_data<T>& init_data, int k, int direction)
        {
            if constexpr (Fields == DataFields::Conservative_Old && std::is_same_v<InsertData, Kokkos::Array<T, 4>>)
            {
                cellsdata.rho_old(k) = data[0];
                cellsdata.Momentum_u_old(k) = data[1];
                cellsdata.Momentum_v_old(k) = data[2];
                cellsdata.E_old(k) = data[3];
                Update_DataField<T, DataFields::DerivedFields>(cellsdata, init_data, k, direction);
            }
            else if constexpr (Fields == DataFields::Conservative_New && std::is_same_v<InsertData, Kokkos::Array<T, 4>>)
            {
                cellsdata.rho_new(k) = data[0];
                cellsdata.Momentum_u_new(k) = data[1];
                cellsdata.Momentum_v_new(k) = data[2];
                cellsdata.E_new(k) = data[3];
            }
            else
            {
                static_assert((Fields == DataFields::Conservative_Old || Fields == DataFields::Conservative_New) && std::is_same_v<InsertData, Kokkos::Array<T, 4>>, "This overload accepts conservative states only");
            }
        }

        template <typename T, DataFields Fields>
        KOKKOS_INLINE_FUNCTION
        void Update_DataField(const CellsData<T>& cellsdata, int k_cell, int direction)
        {
            static_assert(Fields == DataFields::Fluxes, "This overload computes fluxes from the current conservative state");
            auto velocity = extract_DataField<T, DataFields::Velocity>(cellsdata, k_cell);
            auto normal = extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k_cell, direction);
            auto Conservative_state = extract_DataField<T, DataFields::Conservative_Old>(cellsdata, k_cell);
            T Pressure = cellsdata.Pressure(k_cell);
            T Vn = dot(velocity, normal);

            Kokkos::Array<T, 4> Pressure_Energy = { 0, normal[0] * Pressure, normal[1] * Pressure, Pressure * Vn };
            auto flux = Conservative_state * Vn + Pressure_Energy;
            Update_DataField<T, DataFields::Fluxes>(flux, cellsdata, k_cell, direction);
        }

        template <typename ExecSpace, CellCopyFields Fields = CellCopyFields::All, typename T>
        void ApplyPeriodicity(const CellsData<T>& cellsdata)
        {
            Kokkos::parallel_for("Periodicity", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_ghost_layers* (cellsdata.Nj_total_cells)), KOKKOS_LAMBDA(int K_g)
            {
                int i_g = K_g % cellsdata.N_ghost_layers;
                int j_g = K_g / cellsdata.N_ghost_layers;

                const int k_Lcell_destination = Cell_index(i_g, j_g, cellsdata.Ni_total_cells);
                const int k_Lcell_source = k_Lcell_destination + cellsdata.Ni_physical_cells;
                const int k_Rcell_source = k_Lcell_destination + cellsdata.N_ghost_layers;
                const int k_Rcell_destination = k_Rcell_source + cellsdata.Ni_physical_cells;

                copy_cell <Fields>(cellsdata, k_Lcell_source, k_Lcell_destination);
                copy_cell <Fields>(cellsdata, k_Rcell_source, k_Rcell_destination); 
            });
        }
		
        template<typename ExecSpace, typename T>
        void ApplyInitialConditions(const CellsData<T>& cellsdata, const initialization_data<T>& init_data)
        {
            Kokkos::parallel_for("InitialConditions", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_total_cells), KOKKOS_LAMBDA(int k)
            {
                cellsdata.rho_old(k) = init_data.rho_inf;
				cellsdata.velocity_u(k) = init_data.velocity_u_inf;
				cellsdata.velocity_v(k) = init_data.velocity_v_inf;
                cellsdata.Momentum_u_old(k) = init_data.Momentum_u_inf;
                cellsdata.Momentum_v_old(k) = init_data.Momentum_v_inf;
                cellsdata.E_old(k) = init_data.E_inf;
                cellsdata.Pressure(k) = init_data.P_inf;
                cellsdata.SoundSpeed(k) = init_data.c_inf;
                const int i = k % (cellsdata.Ni_total_cells);
                const int j = k / (cellsdata.Ni_total_cells);
                if (i >= cellsdata.N_ghost_layers && i < cellsdata.N_ghost_layers + cellsdata.Ni_physical_cells && j >= cellsdata.N_ghost_layers && j < cellsdata.N_ghost_layers + cellsdata.Nj_physical_cells)
                {
                    const Kokkos::Array<T, 2> Velocity = {init_data.velocity_u_inf, init_data.velocity_v_inf};
                    for (int direction = 0; direction < cellsdata.N_faces; ++direction)
                    {
                        const auto normal = extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k, direction);
                        const T EdgeLength = extract_DataField<T, DataFields::EdgeLength>(cellsdata, k, direction);
                        cellsdata.Lambdas(k, direction) = eval_Lambda<T>(Velocity, normal, init_data.c_inf, EdgeLength);
                    }
                }
            });
        }

        template<typename ExecSpace, typename T, SolverTypes SType = SolverTypes::Euler2D, BoundaryCondition BType = BoundaryCondition::AdiabaticWall>
        void ApplyWallConditions(const CellsData<T>& cellsdata, const initialization_data<T>& init_data, T Wall_temp = T{1})
        {
            Kokkos::parallel_for("WallConditions", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.Ni_physical_cells * cellsdata.N_ghost_layers), KOKKOS_LAMBDA(int k_g)
            {
                int i_g = k_g % cellsdata.Ni_physical_cells;
                int j_g = k_g / cellsdata.Ni_physical_cells;

                const int row_stride = cellsdata.Ni_total_cells;
                const int k_wall = Cell_index(i_g + cellsdata.N_ghost_layers, cellsdata.N_ghost_layers, cellsdata.Ni_total_cells);
                const int k_source = k_wall + (cellsdata.N_ghost_layers - 1 - j_g) * row_stride;
                const int k_destination = k_wall - (cellsdata.N_ghost_layers - j_g) * row_stride;
                const int k_ghost = k_destination;
                
                //eval velocity reflection
                auto ghost_normal = -1 * extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k_wall, 0);
                auto RealCell_velocity = extract_DataField<T, DataFields::Velocity>(cellsdata, k_source);
                auto Ghost_State = extract_DataField<T, DataFields::Conservative_Old>(cellsdata, k_source);
                Kokkos::Array<T, 2> GhostCell_velocity;
                const T wall_temperature = Wall_temp;
                const T gamma = init_data.gamma;
                if constexpr (SType == SolverTypes::Euler2D)
                {
                    GhostCell_velocity = RealCell_velocity - (T)2 * dot(RealCell_velocity, ghost_normal) * ghost_normal;
                }
                if constexpr (SType == SolverTypes::RANS2D)
                {
					
                    if constexpr (BType == BoundaryCondition::AdiabaticWall)
                    {
                        GhostCell_velocity = -1 * RealCell_velocity;
                    }
                    else if constexpr (BType == BoundaryCondition::isothermalWall)
                    {
                        GhostCell_velocity = -1 * RealCell_velocity;
                        const T Ghost_temperature = T{2} * wall_temperature - extract_DataField<T, DataFields::Temperature>(cellsdata, k_source);
                        Ghost_State[0] = cellsdata.Pressure(k_source) / Ghost_temperature;
                        Ghost_State[3] = eval_Energy<T>(Ghost_State[0], GhostCell_velocity[0], GhostCell_velocity[1], cellsdata.Pressure(k_source), gamma);
                    }
                }
                Ghost_State[1] = Ghost_State[0] * GhostCell_velocity[0];
                Ghost_State[2] = Ghost_State[0] * GhostCell_velocity[1];
                const T EdgeLength = extract_DataField<T, DataFields::EdgeLength>(cellsdata, k_wall, 0);
                Update_DataField<T, DataFields::NormalUnitVector>(ghost_normal, cellsdata, k_ghost, 2);
                Update_DataField<T, DataFields::EdgeLength>(EdgeLength, cellsdata, k_ghost, 2);
                // Ghost adjacent to the wall: metrics of the common face seen from the ghost side
                // (both local axes are reversed, so the terms are the opposite of the wall cell's).
                if (k_source == k_wall)
                {
                    for (int component = 0; component < 4; ++component)
                        cellsdata.DJacobien(k_ghost, 2, component) = -cellsdata.DJacobien(k_wall, 0, component);
                }
                Update_DataField<T, DataFields::Conservative_Old>(Ghost_State, cellsdata, init_data, k_ghost, 2);
                Update_DataField<T, DataFields::Fluxes>(cellsdata, k_ghost, 2);
            });
        }
        
        template <typename ExecSpace, typename T>
        void ApplyRiemanInvariantConditions(const CellsData<T>& cellsdata, const initialization_data<T>& init_data)
        {
            Kokkos::parallel_for("RiemannInvariantConditions", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.Ni_physical_cells), KOKKOS_LAMBDA(int i_g)
            {
                int i_cells = i_g + cellsdata.N_ghost_layers;
                int j_Cell = cellsdata.N_ghost_layers + cellsdata.Nj_physical_cells - 1;

                const int row_stride = cellsdata.Ni_total_cells;
                const int k_cell = Cell_index(i_cells, j_Cell, cellsdata.Ni_total_cells);

                Kokkos::Array<T, 4> Ghost_State;
                const auto normal = extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k_cell, 2);
                const auto V_in = extract_DataField<T, DataFields::Velocity>(cellsdata, k_cell);
                const auto V_inf = Kokkos::Array<T, 2>{ init_data.velocity_u_inf, init_data.velocity_v_inf };

                const T c_in = cellsdata.SoundSpeed(k_cell);
                const T c_inf = init_data.c_inf;
                const T rho_in = cellsdata.rho_old(k_cell);
                const T rho_inf = init_data.rho_inf;
                const T gamma = init_data.gamma;
                const T Mach_inf = init_data.Mach_inf;

                T Entropie_k;
                Kokkos::Array<T, 2> V_boundary;

                T Vn_in = dot(V_in, normal);
                T Vn_inf = dot(V_inf, normal);

                T Mach = eval_MachNumber(Vn_in, c_in);

                if (Mach < -1.0)
                {
                    Ghost_State = { init_data.rho_inf, init_data.Momentum_u_inf, init_data.Momentum_v_inf, init_data.E_inf };
                }
                else if (Mach < 1)
                {
                    T J_plus = Vn_in + 2 * c_in / (init_data.gamma - 1.0);
                    T J_minus = Vn_inf - 2 * c_inf / (init_data.gamma - 1.0);
                    T Vbn = (J_plus + J_minus) / T{ 2 };
                    T cb = (gamma - T{ 1 }) / T{ 4 } * (J_plus - J_minus);

                    if (Vbn > T{ 0 })
                    {

                        Entropie_k = c_in * c_in / (gamma * Kokkos::pow(rho_in, gamma - T{ 1 }));
                        V_boundary = V_in + (Vbn - Vn_in) * normal;
                    }
                    else
                    {
                        Entropie_k = c_inf * c_inf / (gamma * Kokkos::pow(rho_inf, gamma - T{ 1 }));
                        V_boundary = V_inf + (Vbn - Vn_inf) * normal;
                    }
                    T rho_boundary = Kokkos::pow(cb * cb / (gamma * Entropie_k), T{ 1 } / (gamma - T{ 1 }));
                    T pressure_boundary = rho_boundary * cb * cb / gamma;

                    Ghost_State = { rho_boundary, rho_boundary * V_boundary[0], rho_boundary * V_boundary[1], eval_Energy<T>(rho_boundary, V_boundary[0], V_boundary[1], pressure_boundary, gamma) };
                }
                else
                {
                    Ghost_State = { cellsdata.rho_old(k_cell), cellsdata.Momentum_u_old(k_cell), cellsdata.Momentum_v_old(k_cell), cellsdata.E_old(k_cell) };
                }
                const auto ghost_normal = T{-1} * normal;
                const T EdgeLength = extract_DataField<T, DataFields::EdgeLength>(cellsdata, k_cell, 2);
                // Ghost adjacent to the farfield: metrics of the common face seen from the ghost side
                // (both local axes are reversed, so the terms are the opposite of the boundary cell's).
                for (int component = 0; component < 4; ++component)
                    cellsdata.DJacobien(k_cell + row_stride, 0, component) = -cellsdata.DJacobien(k_cell, 2, component);
                for (int layer = 1; layer <= cellsdata.N_ghost_layers; layer++)
                {
                    const int k_ghost = k_cell + layer * row_stride;
                    Update_DataField<T, DataFields::NormalUnitVector>(ghost_normal, cellsdata, k_ghost, 0);
                    Update_DataField<T, DataFields::EdgeLength>(EdgeLength, cellsdata, k_ghost, 0);
                    Update_DataField<T, DataFields::Conservative_Old>(Ghost_State, cellsdata, init_data, k_ghost, 0);
                    Update_DataField<T, DataFields::Fluxes>(cellsdata, k_ghost, 0);
                }
            });
        }


        // Complete physical boundaries before extending their ghost rows periodically.
        template <typename ExecSpace, typename T, SolverTypes SType = SolverTypes::Euler2D>
        void ApplyBoundaryConditions(const CellsData<T>& cellsdata, const initialization_data<T>& init_data)
        {
            if (init_data.wall_condition == BoundaryCondition::isothermalWall)
                ApplyWallConditions<ExecSpace, T, SType, BoundaryCondition::isothermalWall>(cellsdata, init_data, init_data.wall_temperature);
            else
                ApplyWallConditions<ExecSpace, T, SType, BoundaryCondition::AdiabaticWall>(cellsdata, init_data, init_data.wall_temperature);
            ApplyRiemanInvariantConditions<ExecSpace>(cellsdata, init_data);
            ApplyPeriodicity<ExecSpace>(cellsdata);
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

        template<typename T>
        KOKKOS_FUNCTION
        Kokkos::Array<T, 4> eval_CellDJacobien(CellsData<T> cellsdata, int k_cell, int direction)
        {
            int row_stride = cellsdata.Ni_total_cells;

            int k_cell_T = Cell_index_sides(k_cell, row_stride, (direction + 1) % 4);
            int k_cell_B = Cell_index_sides(k_cell, row_stride, (direction + 3) % 4);

            int k_side = Cell_index_sides(k_cell, row_stride, direction);
            int k_side_T = Cell_index_sides(k_side, row_stride, (direction + 1) % 4);
            int k_side_B = Cell_index_sides(k_side, row_stride, (direction + 3) % 4);
           
            auto xy_ceta_cell = (extract_DataField<T, DataFields::CellCenter>(cellsdata, k_cell_T) - extract_DataField<T, DataFields::CellCenter>(cellsdata, k_cell_B))/ T{ 2 };
            auto xy_ceta_side = (extract_DataField<T, DataFields::CellCenter>(cellsdata, k_side_T) - extract_DataField<T, DataFields::CellCenter>(cellsdata, k_side_B)) / T{ 2 };
            
            auto xy_cepsilon = extract_DataField<T, DataFields::CellCenter>(cellsdata, k_side) - extract_DataField<T, DataFields::CellCenter>(cellsdata, k_cell);
            auto xy_ceta = mean<T>(xy_ceta_cell, xy_ceta_side);

            T J_f = det2D(xy_cepsilon, xy_ceta);
            return Kokkos::Array<T, 4>{xy_cepsilon[0]/J_f, xy_cepsilon[1]/J_f, xy_ceta[0]/J_f, xy_ceta[1]/J_f};
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
        Kokkos::Array<T, 3> eval_CellCenterAreas(int i, int j, const MeshData_Structured<T>& meshdata)
        {
            int Ni = meshdata.Ni;

			//Positions have {x,y} coordinates
			Kokkos::Array <T, 2> P_1 = { meshdata.X(Node_index(i, j, Ni)), meshdata.Y(Node_index(i, j, Ni)) };
			Kokkos::Array <T, 2> P_2 = { meshdata.X(Node_index(i + 1, j, Ni)), meshdata.Y(Node_index(i + 1, j, Ni)) };
			Kokkos::Array <T, 2> P_3 = { meshdata.X(Node_index(i + 1, j + 1, Ni)), meshdata.Y(Node_index(i + 1, j + 1, Ni)) };
			Kokkos::Array <T, 2> P_4 = { meshdata.X(Node_index(i, j + 1, Ni)), meshdata.Y(Node_index(i, j + 1, Ni)) };
            

			auto S_1 = P_4 - P_1;
			auto S_2 = P_3 - P_1;
			auto S_3 = P_2 - P_1;


			T Cell_Area1 = T{ 0.5 } * det2D(S_2, S_1);
			T Cell_Area2 = T{ 0.5 } * det2D(S_3, S_2);
			T Cell_Area = Cell_Area1 + Cell_Area2;

			auto P_center1 = (P_1 + P_3 + P_4) / T{ 3 };
			auto P_center2 = (P_1 + P_2 + P_3) / T{ 3 };
			auto P_center = (P_center1 * Cell_Area1 + P_center2 * Cell_Area2) / Cell_Area;

			Kokkos::Array<T, 3> CellCenterAreas = { P_center[0], P_center[1], Cell_Area };
            return CellCenterAreas;
        }

        template <typename ExecSpace, typename T>
        void eval_GeometricalData(const CellsData<T>& cellsdata, const MeshData_Structured<T>& meshdata)
        {
			
            Kokkos::parallel_for("eval_inner_geometrical_data", Kokkos::RangePolicy<ExecSpace>(0, meshdata.N_Cells), KOKKOS_LAMBDA(int k_mesh)
            {
                int i = k_mesh % cellsdata.Ni_physical_cells;
                int j = k_mesh / cellsdata.Ni_physical_cells;

                int i_cell = i + cellsdata.N_ghost_layers;
                int j_cell = j + cellsdata.N_ghost_layers;
                int k_cell = Cell_index(i_cell, j_cell, cellsdata.Ni_total_cells);
				auto CellCenterAreas = eval_CellCenterAreas(i, j, meshdata);
                cellsdata.Cell_Areas(k_cell) = CellCenterAreas[2];
				cellsdata.Cell_Centers(k_cell, 0) = CellCenterAreas[0];
				cellsdata.Cell_Centers(k_cell, 1) = CellCenterAreas[1];

                for (int direction = 0; direction < cellsdata.N_faces; direction++)
                {
                    const auto normals = eval_CellNormal(i, j, direction, meshdata);
                    cellsdata.Normals(k_cell, direction, 0) = normals[0];
                    cellsdata.Normals(k_cell, direction, 1) = normals[1];
                    cellsdata.Normals(k_cell, direction, 2) = normals[2];                    
                }

            });
            Kokkos::parallel_for("eval_ghost_wall_center_data", Kokkos::RangePolicy<ExecSpace>(0, (cellsdata.Ni_physical_cells*cellsdata.N_ghost_layers)), KOKKOS_LAMBDA(int k_mesh)
            {
				int row_stride = cellsdata.Ni_total_cells;

                int i_mesh = k_mesh % cellsdata.Ni_physical_cells ;
                int j_mesh = k_mesh / cellsdata.Ni_physical_cells ;

                int i_cell = i_mesh + cellsdata.N_ghost_layers;
                int j_cell = j_mesh + cellsdata.N_ghost_layers;

                int k_wall_node = Node_index(i_mesh, 0, cellsdata.Ni_nodes);

                int k_wall = Cell_index(i_cell, cellsdata.N_ghost_layers, cellsdata.Ni_total_cells);
                int k_cell = k_wall + j_mesh * row_stride;
                int k_ghost = k_wall - (j_mesh + 1) * row_stride;


                
                Kokkos::Array<T, 2> P_reference = {meshdata.X(k_wall_node), meshdata.Y(k_wall_node) };
				Kokkos::Array<T, 2> P_F1 = { meshdata.X(k_wall_node + 1), meshdata.Y(k_wall_node + 1)};
				auto P_Cf_relative = (P_F1 - P_reference)/ T{ 2 };
                auto P_Ccell_relative = (extract_DataField<T, DataFields::CellCenter>(cellsdata, k_cell)) - P_reference;
				auto P_Cghost_relative = 2 * P_Cf_relative - P_Ccell_relative;
				auto P_Cghost = P_Cghost_relative + P_reference;

                cellsdata.Cell_Centers(k_ghost, 0) = P_Cghost[0];
                cellsdata.Cell_Centers(k_ghost, 1) = P_Cghost[1];
            });
            Kokkos::parallel_for("eval_ghost_Farfield_center_data", Kokkos::RangePolicy<ExecSpace>(0, (cellsdata.Ni_physical_cells* cellsdata.N_ghost_layers)), KOKKOS_LAMBDA(int k_mesh)
            {
                int row_stride = cellsdata.Ni_total_cells;
                int i_mesh = k_mesh % cellsdata.Ni_physical_cells;
                int j_mesh = k_mesh / cellsdata.Ni_physical_cells;


                int i_cell = i_mesh + cellsdata.N_ghost_layers;
				int j_farfield = cellsdata.N_ghost_layers + cellsdata.Nj_physical_cells;
                int k_Farfield = Cell_index(i_cell, j_farfield, cellsdata.Ni_total_cells);
				int k_Farfirld_node = Node_index(i_mesh, cellsdata.Nj_nodes - 1, cellsdata.Ni_nodes);
                int k_cell = k_Farfield - (cellsdata.N_ghost_layers - j_mesh)  * row_stride;
                int k_ghost = k_Farfield + (cellsdata.N_ghost_layers- 1 - j_mesh) * row_stride;
                



                Kokkos::Array<T, 2> P_reference = { meshdata.X(k_Farfirld_node), meshdata.Y(k_Farfirld_node) };
                Kokkos::Array<T, 2> P_F1 = { meshdata.X(k_Farfirld_node + 1), meshdata.Y(k_Farfirld_node + 1) };
                auto P_Cf_relative = (P_F1 - P_reference) / T{ 2 };
                auto P_Ccell_relative = (extract_DataField<T, DataFields::CellCenter>(cellsdata, k_cell)) - P_reference;
                auto P_Cghost_relative = 2 * P_Cf_relative - P_Ccell_relative;
                auto P_Cghost = P_Cghost_relative + P_reference;

                cellsdata.Cell_Centers(k_ghost, 0) = P_Cghost[0];
                cellsdata.Cell_Centers(k_ghost, 1) = P_Cghost[1];
            });

            ApplyPeriodicity<ExecSpace, CellCopyFields::CellCenter>(cellsdata);
            ExecSpace().fence("Cell centers ready for metrics");
            Kokkos::parallel_for("eval_inner_DJacobien", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_physical_cells), KOKKOS_LAMBDA(int k_mesh)
            {
                const int i_cell = k_mesh % cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
                const int j_cell = k_mesh / cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
                const int k_cell = Cell_index(i_cell, j_cell, cellsdata.Ni_total_cells);
                for (int direction = 0; direction < cellsdata.N_faces; ++direction)
                {
                    const auto DJacobien = eval_CellDJacobien(cellsdata, k_cell, direction);
                    for (int component = 0; component < 4; ++component)
                        cellsdata.DJacobien(k_cell, direction, component) = DJacobien[component];
                }
            });
        }


        template <typename T, DataFields Field>
        KOKKOS_INLINE_FUNCTION
        auto eval_gradient(const CellsData<T>& cellsdata, int k_cell, int direction)
        {
			int row_stride = cellsdata.Ni_total_cells;
            
			int k_cell_T = Cell_index_sides(k_cell, row_stride, (direction + 1) % 4);
			int k_cell_B = Cell_index_sides(k_cell, row_stride, (direction + 3) % 4);

            int k_side = Cell_index_sides(k_cell, row_stride, direction);
            int k_side_T = Cell_index_sides(k_side, row_stride, (direction + 1) % 4);
            int k_side_B = Cell_index_sides(k_side, row_stride, (direction + 3) % 4);
	
            
			auto DA_eta_cell = (extract_DataField<T, Field>(cellsdata, k_cell_T) - extract_DataField<T, Field>(cellsdata, k_cell_B)) / T{ 2 };
			auto DA_eta_side = (extract_DataField<T, Field>(cellsdata, k_side_T) - extract_DataField<T, Field>(cellsdata, k_side_B)) / T{ 2 };
            
            
            auto DA_epsilon = (extract_DataField<T, Field>(cellsdata, k_side) - extract_DataField<T, Field>(cellsdata, k_cell));
            auto DA_eta = mean<T>(DA_eta_cell, DA_eta_side);

            auto DAx = cellsdata.DJacobien(k_cell, direction, 3) * DA_epsilon - cellsdata.DJacobien(k_cell, direction, 1)*DA_eta;
            auto DAy = cellsdata.DJacobien(k_cell, direction, 0) * DA_eta - cellsdata.DJacobien(k_cell, direction, 2)*DA_epsilon;


            if constexpr (Field == DataFields::Velocity)
            {
                return Kokkos::Array<T, 4>{DAx[0], DAx[1], DAy[0], DAy[1]};
            }
			else if constexpr (Field == DataFields::Temperature)
            {
                return Kokkos::Array<T, 2>{DAx, DAy};
            }
            else
            {
				static_assert(Field == DataFields::Velocity || Field == DataFields::Temperature, "Field must be either Velocity or Temperature");
            }
        }

        template <typename T>
        KOKKOS_INLINE_FUNCTION
        Kokkos::Array<T, 4> eval_viscuousDissipation(const CellsData<T>& cellsdata, const initialization_data<T>& data, int k_cell, int direction)
        {

            int row_stride = cellsdata.Ni_total_cells;
			int k_side = Cell_index_sides(k_cell, row_stride, direction);
			Kokkos::Array<T, 4> ViscuousDissipation;
			auto velocity = extract_DataField<T, DataFields::Velocity>(cellsdata, k_cell);
			auto velocity_side = extract_DataField<T, DataFields::Velocity>(cellsdata, k_side);
            auto normal = extract_DataField<T, DataFields::NormalUnitVector>(cellsdata, k_cell, direction);
			
            
            auto dV = eval_gradient<T, DataFields::Velocity>(cellsdata, k_cell, direction);

            
            auto tau = Tau(data.mu_inf, dV);
			auto t_x = normal[0] * tau[0] + normal[1] * tau[2];
			auto t_y = normal[0] * tau[2] + normal[1] * tau[1];
			auto V_mean = mean(velocity, velocity_side);
			auto dtheta = eval_gradient<T, DataFields::Temperature>(cellsdata, k_cell, direction);
            auto q = -data.k_inf * dot(dtheta, normal);

            ViscuousDissipation = {0, t_x, t_y, V_mean[0]*t_x + V_mean[1]*t_y - q};

            return ViscuousDissipation;
        }

        template <typename T, SolverTypes SType = SolverTypes::Euler2D>
        KOKKOS_INLINE_FUNCTION
        Kokkos::Array<T, 4> eval_Fluxes(const CellsData<T>& cellsdata, const initialization_data<T>& init_data, int k_cell, int direction)
        {

            int row_stride = cellsdata.Ni_total_cells;
            Kokkos::Array<T, 4> Fluxes_2;
                
            int k_side = Cell_index_sides(k_cell, row_stride, direction);
            auto Fluxes_k_cell = extract_DataField<T, DataFields::Fluxes>(cellsdata, k_cell, direction);
            auto Fluxes_k_side = extract_DataField<T, DataFields::Fluxes>(cellsdata, k_side, (direction + 2) % 4);
            if constexpr (SType == SolverTypes::RANS2D)
            {
				auto ViscuousDissipation_cell = eval_viscuousDissipation<T>(cellsdata, init_data, k_cell, direction);
                auto ViscuousDissipation_side = eval_viscuousDissipation<T>(cellsdata, init_data, k_side, (direction + 2) % 4);
				Fluxes_k_cell = Fluxes_k_cell - ViscuousDissipation_cell;
				Fluxes_k_side = Fluxes_k_side - ViscuousDissipation_side;
            }

			Fluxes_2 = mean<T>(Fluxes_k_cell, -1 * Fluxes_k_side);
            return Fluxes_2;
        }

        template <typename T>
        KOKKOS_INLINE_FUNCTION
        Kokkos::Array<T, 4> eval_ArtDissipation(const CellsData<T>& cellsdata, int k_cell, int direction)
        {
			int row_stride = cellsdata.Ni_total_cells;
            
            const T k_2 = (T)0.6;
            const T k_4 = (T)1/(T)64;
            
            T Lambdas_2 = mean(cellsdata.Lambdas(k_cell, direction), cellsdata.Lambdas(Cell_index_sides(k_cell, row_stride, direction), (direction + 2)%4));
			Kokkos::Array <T, 4> Pressure_captors;

            auto k_dissipations = [](int direction, int row_stride, int k_cell)
                {
                    switch (direction)
                    {
                    case 0: return Kokkos::Array<int, 6>{k_cell + 2*row_stride, k_cell + row_stride, k_cell, k_cell - row_stride, k_cell - 2 * row_stride, k_cell - 3 * row_stride};
                          break;
                    case 1: return Kokkos::Array<int, 6>{k_cell - 2,k_cell - 1, k_cell, k_cell + 1, k_cell + 2, k_cell + 3};
                          break;
                    case 2: return Kokkos::Array<int, 6>{k_cell - 2*row_stride,k_cell - row_stride, k_cell, k_cell + row_stride, k_cell + 2 * row_stride, k_cell + 3 * row_stride};
                          break;
                    case 3: return Kokkos::Array<int, 6>{k_cell + 2,k_cell + 1, k_cell, k_cell - 1, k_cell - 2, k_cell - 3};
                          break;
                    }
                };
            Kokkos::Array<int, 6> k_Dissipation = k_dissipations(direction, row_stride, k_cell);
            const auto& Pressure = cellsdata.Pressure;
            for (int Captor_index = 0; Captor_index < 4; Captor_index++)
            {
                int k_dissipation_index = Captor_index + 1;
                T denominator = Pressure(k_Dissipation[k_dissipation_index + 1]) + 2 * Pressure(k_Dissipation[k_dissipation_index]) + Pressure(k_Dissipation[k_dissipation_index - 1]);
                T numerator = abs(Pressure(k_Dissipation[k_dissipation_index + 1]) - 2 * Pressure(k_Dissipation[k_dissipation_index]) + Pressure(k_Dissipation[k_dissipation_index - 1]));
                Pressure_captors[Captor_index] = numerator / denominator;
            }
            
            T epsilon_2 = k_2 * max(Pressure_captors);

            T epsilon_4 = max(Kokkos::Array<T,2> {0 , k_4 - epsilon_2});

            Kokkos::Array < Kokkos::Array<T, 4>, 4> W;
            for (int i = 0; i < 4; i++)
            {
                W[i] = extract_DataField<T, DataFields::Conservative_Old>(cellsdata, k_Dissipation[i + 1]);
            }
            Kokkos::Array<T, 4> Dissipation = Lambdas_2 * (epsilon_2* (W[2]-W[1]) - epsilon_4* (W[3]-3*W[2] + 3*W[1]-W[0]));
            return Dissipation;
        }

        template <typename T, SolverTypes SType = SolverTypes::Euler2D>
        KOKKOS_INLINE_FUNCTION
        Kokkos::pair<Kokkos::Array<T, 4>, Kokkos::Array<T, 4>> eval_All_FluxesDissipation(const CellsData<T>& cellsdata, const initialization_data<T>& init_data, int k_cell, bool evaluate_dissipation = true)
        {
            Kokkos::Array<T, 4> Fluxes = {0, 0, 0, 0};
            Kokkos::Array<T, 4> Dissipation = {0, 0, 0, 0};
            for (int direction = 0; direction < cellsdata.N_faces; ++direction)
            {
                const T edge_length = extract_DataField<T, DataFields::EdgeLength>(cellsdata, k_cell, direction);
                Fluxes += edge_length * eval_Fluxes<T, SType>(cellsdata, init_data, k_cell, direction);
                if (evaluate_dissipation)
                    Dissipation += eval_ArtDissipation<T>(cellsdata, k_cell, direction);
            }
            return Kokkos::make_pair(Fluxes, Dissipation);
        }
     
        template <typename T, SolverTypes SType = SolverTypes::Euler2D>
        KOKKOS_INLINE_FUNCTION
        T eval_TimeStep(const CellsData<T>& cellsdata, const initialization_data<T>& init_data, int k_cell, T CFL)
        {
            T Lambda = T{0};
            for (int direction = 0; direction < cellsdata.N_faces; ++direction)
            {
                Lambda += T{0.5} * cellsdata.Lambdas(k_cell, direction);
                if constexpr (SType == SolverTypes::RANS2D)
                {
                    const T edge_length = extract_DataField<T, DataFields::EdgeLength>(cellsdata, k_cell, direction);
                    Lambda += (init_data.mu_inf / cellsdata.rho_old(k_cell)) * (T{4}/T{3} + init_data.gamma/init_data.Pr_inf) * edge_length*edge_length / cellsdata.Cell_Areas(k_cell);
                }
            }
            return CFL * cellsdata.Cell_Areas(k_cell) / Lambda;
        }

        #pragma region Runge-Kutta shared steps
        // Storage index (ghost layers included) of the k_mesh-th physical cell.
        template <typename T>
        KOKKOS_INLINE_FUNCTION
        int Physical_Cell_index(const CellsData<T>& cellsdata, int k_mesh)
        {
            const int i_cell = k_mesh % cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
            const int j_cell = k_mesh / cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
            return Cell_index(i_cell, j_cell, cellsdata.Ni_total_cells);
        }

        template <typename T>
        void RungeKutta_Validate(const CellsData<T>& cellsdata, T CFL, T max_time_step, bool physical_time_step)
        {
            if (cellsdata.N_faces != 4 || cellsdata.N_ghost_layers < 3)
                throw std::invalid_argument("Runge-Kutta requires four faces and at least three ghost layers");
            if (cellsdata.Ni_physical_cells < cellsdata.N_ghost_layers || cellsdata.Nj_physical_cells < cellsdata.N_ghost_layers)
                throw std::invalid_argument("Runge-Kutta boundary stencils require enough physical cells for the ghost layers");
            if (!(CFL > T{0}))
                throw std::invalid_argument("Runge-Kutta requires a positive CFL");
            if (physical_time_step && (!(max_time_step > T{0}) || !std::isfinite(max_time_step)))
                throw std::invalid_argument("Runge-Kutta requires a positive finite physical step limit");
        }

        // Derived fields and fluxes of the physical cells, then the ghost cells; optionally the local time steps.
        template <typename ExecSpace, typename T, SolverTypes SType>
        void RungeKutta_PrepareFields(const CellsData<T>& cellsdata, const initialization_data<T>& init_data, T CFL, bool evaluate_time_step)
        {
            Kokkos::parallel_for("RungeKutta_PrepareFields", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_physical_cells), KOKKOS_LAMBDA(int k_mesh)
            {
                const int k_cell = Physical_Cell_index(cellsdata, k_mesh);
                Update_DataField<T, DataFields::DerivedFields>(cellsdata, init_data, k_cell);
                for (int direction = 0; direction < cellsdata.N_faces; ++direction)
                    Update_DataField<T, DataFields::Fluxes>(cellsdata, k_cell, direction);
                if (evaluate_time_step)
                    cellsdata.Time_Step(k_cell) = eval_TimeStep<T, SType>(cellsdata, init_data, k_cell, CFL);
            });
            ExecSpace().fence("DerivedFields are ready, preparing ghosts");
            ApplyBoundaryConditions<ExecSpace, T, SType>(cellsdata, init_data);
            ExecSpace().fence("GhostCells are ready, RungeKutta is Ready");
        }

        // Global physical time step: smallest local step, capped by max_time_step.
        template <typename ExecSpace, typename T>
        T RungeKutta_MinTimeStep(const CellsData<T>& cellsdata, T max_time_step)
        {
            T min_time_step = max_time_step;
            Kokkos::parallel_reduce("RungeKutta_MinTimeStep", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_physical_cells), KOKKOS_LAMBDA(int k_mesh, T& local_min_time_step)
            {
                local_min_time_step = Kokkos::min<T>(local_min_time_step, cellsdata.Time_Step(Physical_Cell_index(cellsdata, k_mesh)));
            }, Kokkos::Min<T>(min_time_step));
            return Kokkos::min(min_time_step, max_time_step);
        }

        // The stage result becomes the current state.
        template <typename T>
        void RungeKutta_SwapStates(CellsData<T>& cellsdata)
        {
            std::swap(cellsdata.rho_new, cellsdata.rho_old);
            std::swap(cellsdata.Momentum_u_new, cellsdata.Momentum_u_old);
            std::swap(cellsdata.Momentum_v_new, cellsdata.Momentum_v_old);
            std::swap(cellsdata.E_new, cellsdata.E_old);
        }
        #pragma endregion

        // Jameson hybrid five-stage scheme: dissipation refreshed at stages 1, 3 and 5.
        template <typename ExecSpace, typename T, SolverTypes SType = SolverTypes::Euler2D, TimeEvolutionScheme TEScheme = TimeEvolutionScheme::Steady>
        T RungeKutta5(CellsData<T>& cellsdata, const initialization_data<T>& init_data, T CFL = T{1}, T max_time_step = std::numeric_limits<T>::max())
        {
            RungeKutta_Validate(cellsdata, CFL, max_time_step, TEScheme == TimeEvolutionScheme::UnSteady);

            Kokkos::deep_copy(cellsdata.rho_0, cellsdata.rho_old);
            Kokkos::deep_copy(cellsdata.Momentum_u_0, cellsdata.Momentum_u_old);
            Kokkos::deep_copy(cellsdata.Momentum_v_0, cellsdata.Momentum_v_old);
            Kokkos::deep_copy(cellsdata.Energy_0, cellsdata.E_old);

            const Kokkos::Array<T, 5> Alpha = {T{1}/T{4}, T{1}/T{6}, T{3}/T{8}, T{1}/T{2}, T{1}};
            const Kokkos::Array<T, 5> Beta = {T{1}, T{0}, T(0.56), T{0}, T(0.44)};

            T min_time_step = T{0};

            for (int KuntaStep = 0; KuntaStep < 5; ++KuntaStep)
            {
                RungeKutta_PrepareFields<ExecSpace, T, SType>(cellsdata, init_data, CFL, KuntaStep == 0);
                if constexpr (TEScheme == TimeEvolutionScheme::UnSteady)
                    if (KuntaStep == 0)
                        min_time_step = RungeKutta_MinTimeStep<ExecSpace, T>(cellsdata, max_time_step);

                const bool evaluate_dissipation = KuntaStep == 0 || KuntaStep == 2 || KuntaStep == 4;

                Kokkos::parallel_for("RungeKutta5_UpdateState", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_physical_cells), KOKKOS_LAMBDA(int k_mesh)
                {
                    const int k_cell = Physical_Cell_index(cellsdata, k_mesh);
                    auto Fluxes_Dissipations = eval_All_FluxesDissipation<T, SType>(cellsdata, init_data, k_cell, evaluate_dissipation);
                    auto Fluxes = Fluxes_Dissipations.first;
                    Kokkos::Array<T, 4> Dissipations;

                    if (KuntaStep == 0)
                    {
                        Dissipations = Beta[KuntaStep] * Fluxes_Dissipations.second;
                    }
                    else if (evaluate_dissipation)
                    {
                        auto Previous_Dissipation = extract_DataField<T, DataFields::Dissipation_History>(cellsdata, k_cell);
                        Dissipations = Beta[KuntaStep] * Fluxes_Dissipations.second + (1 - Beta[KuntaStep]) * Previous_Dissipation;
                    }
                    else
                    {
                        Dissipations = extract_DataField<T, DataFields::Dissipation_History>(cellsdata, k_cell);
                    }
                    Update_DataField<T, DataFields::Dissipation_History>(Dissipations, cellsdata, k_cell);

                    auto Residual = Fluxes - Dissipations;
                    Kokkos::Array<T, 4> W_0 = { cellsdata.rho_0(k_cell), cellsdata.Momentum_u_0(k_cell),cellsdata.Momentum_v_0(k_cell), cellsdata.Energy_0(k_cell) };
                    const T time_step = TEScheme == TimeEvolutionScheme::UnSteady ? min_time_step : cellsdata.Time_Step(k_cell);
                    const auto W_new = W_0 - Alpha[KuntaStep] * time_step * Residual / cellsdata.Cell_Areas(k_cell);
                    Update_DataField<T, DataFields::Conservative_New>(W_new, cellsdata, init_data, k_cell);
                });

                RungeKutta_SwapStates(cellsdata);
            }
            // Return with derived fields and ghost cells matching the final state.
            RungeKutta_PrepareFields<ExecSpace, T, SType>(cellsdata, init_data, CFL, false);
            return min_time_step; // Physical time increment; zero for local steady marching.
        }

        // Carpenter-Kennedy fourth-order, five-stage, low-storage (2N) scheme:
        //   G = A[stage] * G + dt * L(W),   W = W + B[stage] * G,   with L(W) = -(Fluxes - Dissipation) / Area.
        // Only the current state and the accumulator G (kept in D_history) are stored.
        template <typename ExecSpace, typename T, SolverTypes SType = SolverTypes::Euler2D, TimeEvolutionScheme TEScheme = TimeEvolutionScheme::UnSteady>
        T RungeKuttaCarpenterKennedy5(CellsData<T>& cellsdata, const initialization_data<T>& init_data, T CFL = T{1}, T max_time_step = std::numeric_limits<T>::max())
        {
            static_assert(TEScheme == TimeEvolutionScheme::UnSteady, "Carpenter-Kennedy is only available in UnSteady mode");
            RungeKutta_Validate(cellsdata, CFL, max_time_step, true);

            const Kokkos::Array<T, 5> A = { T{ 0 }, -T(567301805773LL) / T(1357537059087LL), -T(2404267990393LL) / T(2016746695238LL), -T(3550918686646LL) / T(2091501179385LL), -T(1275806237668LL) / T(842570457699LL) };
            const Kokkos::Array<T, 5> B = { T(1432997174477LL) / T(9575080441755LL), T(5161836677717LL) / T(13612068292357LL), T(1720146321549LL) / T(2090206949498LL), T(3134564353537LL) / T(4481467310338LL), T(2277821191437LL) / T(14882151754819LL) };

            T time_step = T{0};

            for (int KuntaStep = 0; KuntaStep < 5; ++KuntaStep)
            {
                RungeKutta_PrepareFields<ExecSpace, T, SType>(cellsdata, init_data, CFL, KuntaStep == 0);
                if (KuntaStep == 0)
                    time_step = RungeKutta_MinTimeStep<ExecSpace, T>(cellsdata, max_time_step);

                Kokkos::parallel_for("RungeKuttaCarpenterKennedy5_UpdateState", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_physical_cells), KOKKOS_LAMBDA(int k_mesh)
                {
                    const int k_cell = Physical_Cell_index(cellsdata, k_mesh);
                    const auto Fluxes_Dissipations = eval_All_FluxesDissipation<T, SType>(cellsdata, init_data, k_cell);
                    const auto Residual = Fluxes_Dissipations.first - Fluxes_Dissipations.second;
                    // The first stage starts a new accumulator: nothing is carried over from the previous time step.
                    Kokkos::Array<T, 4> G_accumulator = T{-1} * time_step * Residual / cellsdata.Cell_Areas(k_cell);
                    if (KuntaStep > 0)
                        G_accumulator = A[KuntaStep] * extract_DataField<T, DataFields::Dissipation_History>(cellsdata, k_cell) + G_accumulator;
                    const auto W_old = extract_DataField<T, DataFields::Conservative_Old>(cellsdata, k_cell);
                    const auto W_new = W_old + B[KuntaStep] * G_accumulator;
                    Update_DataField<T, DataFields::Dissipation_History>(G_accumulator, cellsdata, k_cell);
                    Update_DataField<T, DataFields::Conservative_New>(W_new, cellsdata, init_data, k_cell);
                });

                RungeKutta_SwapStates(cellsdata);
            }
            // Return with derived fields and ghost cells matching the final state.
            RungeKutta_PrepareFields<ExecSpace, T, SType>(cellsdata, init_data, CFL, false);
            return time_step; // Physical time increment.
        }
}
