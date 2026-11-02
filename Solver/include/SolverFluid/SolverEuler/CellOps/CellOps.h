#pragma once

#include "../CellsData.h"
#include "CellCopy.h"

namespace Solver::CellOps
{
    KOKKOS_INLINE_FUNCTION
    int Node_index(int i, int j, int Ni) {return i + Ni*j;}

    KOKKOS_INLINE_FUNCTION
    int Cell_index(int i, int j, int Ni) {return i + (Ni-1)*j;}

    /**
    * Evaluates the geometrical data for a given cell.
    * @param i The cell index in the i-direction.
    * @param j The cell index in the j-direction.
    * @param direction The direction of the normal vector.
    * @param meshdata The mesh data.
    */
    KOKKOS_INLINE_FUNCTION
    Kokkos::pair<float, float> eval_Cell_Normal(int i, int j, int direction, const MeshData_Euler& meshdata)
    {
        int Ni = meshdata.Ni;
        float delta_x = 0.0f, delta_y = 0.0f;

        switch (direction)
        {
            case 0: // Bottom
                delta_x = meshdata.X(Node_index(i + 1,j, Ni)) - meshdata.X(Node_index(i,j, Ni));
                delta_y = meshdata.Y(Node_index(i + 1,j, Ni)) - meshdata.Y(Node_index(i,j, Ni));

                break;

            case 1: // Right
                delta_x = meshdata.X(Node_index(i+1,j + 1, Ni)) - meshdata.X(Node_index(i+1,j, Ni));
                delta_y = meshdata.Y(Node_index(i+1,j + 1, Ni)) - meshdata.Y(Node_index(i+1,j, Ni));

                break;
            case 2: // Top
                delta_x = meshdata.X(Node_index(i,j + 1, Ni)) - meshdata.X(Node_index(i + 1,j + 1, Ni));
                delta_y = meshdata.Y(Node_index(i,j + 1, Ni)) - meshdata.Y(Node_index(i + 1,j + 1, Ni));

                break;
            case 3: // Left
                delta_x = meshdata.X(Node_index(i,j, Ni)) - meshdata.X(Node_index(i,j + 1, Ni));
                delta_y = meshdata.Y(Node_index(i,j, Ni)) - meshdata.Y(Node_index(i,j + 1, Ni));

                break;
        }
        return {delta_y, - delta_x};
    }

    /**
    * Evaluates the area of a given cell.
    * @param i The cell index in the i-direction.
    * @param j The cell index in the j-direction.
    * @param meshdata The mesh data.
    * @return The area of the cell.
    */
    KOKKOS_INLINE_FUNCTION
    float eval_Cell_Area(int i, int j, const MeshData_Euler& meshdata)
    {
        int Ni = meshdata.Ni;
        float delta_x1 = meshdata.X(Node_index(i + 1, j+ 1, Ni)) - meshdata.X(Node_index(i, j, Ni));
        float delta_y1 = meshdata.Y(Node_index(i + 1, j + 1, Ni)) - meshdata.Y(Node_index(i, j, Ni));
            
        float delta_x2 = meshdata.X(Node_index(i , j + 1, Ni)) - meshdata.X(Node_index(i + 1, j, Ni));
        float delta_y2 = meshdata.Y(Node_index(i, j + 1, Ni)) - meshdata.Y(Node_index(i + 1, j, Ni));

        return 0.5f * ((delta_x1 * delta_y2 - delta_x2 * delta_y1));
    }

    template <typename ExecSpace = Kokkos::DefaultExecutionSpace>
    KOKKOS_FUNCTION
	void eval_geometrical_data(CellsData cellsdata, const MeshData_Euler& meshdata)
	{
        Kokkos::parallel_for("eval_geometrical_data", Kokkos::RangePolicy<ExecSpace>(0, meshdata.N_Cells), KOKKOS_LAMBDA(int k_origin)
        {
            int i = k_origin % (meshdata.Ni - 1);
            int j = k_origin / (meshdata.Ni - 1);

			int i_cell = i + cellsdata.N_gCell_layers;
			int j_cell = j + cellsdata.N_gCell_layers;
            int k_cell = Cell_index(i_cell, j_cell, cellsdata.Nitotal);
            cellsdata.Cell_Areas(k_cell) = eval_Cell_Area(i,j, meshdata);


            for (int direction = 0; direction < meshdata.N_edges; direction++)
            {
                cellsdata.Normals(k_cell, direction, 0) = eval_Cell_Normal(i, j, direction, meshdata).first;
                cellsdata.Normals(k_cell, direction, 1) = eval_Cell_Normal(i, j, direction, meshdata).second;
            }
        });
        Kokkos::fence();
	}

    template <typename ExecSpace = Kokkos::DefaultExecutionSpace,
                CellCopyFields Fields = CellCopyFields::All>
    KOKKOS_FUNCTION
    void ApplyPeriodicity(CellsData cellsdata, const MeshData_Euler& meshdata)
    {
        Kokkos::parallel_for("Periodicity", Kokkos::RangePolicy<ExecSpace>(0, meshdata.Nj-1), KOKKOS_LAMBDA(int j)
        {

			int j_cell = j + cellsdata.N_gCell_layers;
            for (int i_cell_destination = 0; i_cell_destination < cellsdata.N_gCell_layers; i_cell_destination++)
            {
                int i_source = (meshdata.Ni - 1) + i_cell_destination;
				int k_cell_source = Cell_index(i_source, j_cell, cellsdata.Nitotal);
                int k_cell_destination = Cell_index(i_cell_destination, j_cell, cellsdata.Nitotal);
                copy_cell <Fields>(cellsdata, k_cell_source, k_cell_destination);
            }
            for (int i_cell_destination = (cellsdata.Nitotal-1)-cellsdata.N_gCell_layers; i_cell_destination < cellsdata.Nitotal-1; i_cell_destination++)
            {
                int i_source = i_cell_destination - (meshdata.Ni - 1);
                int k_cell_source = Cell_index(i_source, j_cell, cellsdata.Nitotal);
                int k_cell_destination = Cell_index(i_cell_destination, j_cell, cellsdata.Nitotal);
                copy_cell <Fields>(cellsdata, k_cell_source, k_cell_destination);
            }
        });
        Kokkos::fence();
    }
}
