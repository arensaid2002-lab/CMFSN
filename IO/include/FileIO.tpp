#include <sstream>
#pragma once

#include <Kokkos_Core.hpp>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <locale>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <MeshHandler.h>

using namespace std;
using namespace Mesh;

namespace IO
{
    inline std::filesystem::path Create_RunDirectory(const std::filesystem::path& root)
    {
        std::filesystem::create_directories(root);
        for (unsigned long long run = 1; ; ++run)
        {
            std::ostringstream name;
            name << "run_" << std::setfill('0') << std::setw(4) << run;
            const auto directory = root / name.str();
            std::error_code error;
            if (std::filesystem::create_directory(directory, error))
                return directory;
            if (error && error != std::errc::file_exists)
                throw std::filesystem::filesystem_error("Could not create simulation output directory", directory, error);
        }
    }

    template <typename T, MeshFileFormat Format>
    void Import_Mesh(const string& path, MeshData_Structured<T>& MeshData)
    {
        try
        {
            ifstream Meshfile(path);
            if (Meshfile.is_open())
            {
                MeshData_Structured<T> imported;
                if (!(Meshfile >> imported.Ni >> imported.Nj) || imported.Ni < 2 || imported.Nj < 2)
                    throw runtime_error("Invalid mesh dimensions");
                if (imported.Ni > numeric_limits<int>::max() / imported.Nj)
                    throw runtime_error("Mesh dimensions exceed the supported index range");
                imported.N_nodes = imported.Ni * imported.Nj;
                imported.N_Cells = (imported.Ni-1) * (imported.Nj-1);
 
                Kokkos::View<T*, Kokkos::LayoutLeft > X("X", imported.N_nodes);
                Kokkos::View<T*, Kokkos::LayoutLeft > Y("Y", imported.N_nodes);

                auto X_host = Kokkos::create_mirror_view(X);
                auto Y_host = Kokkos::create_mirror_view(Y);

                for(int i = 0; i < imported.N_nodes; i++)
                {
                    if(!(Meshfile >> X_host(i)))
                    {
                        throw runtime_error("X Coordinate reading failed");
                    }
                }

                for(int i = 0; i < imported.N_nodes; i++)
                {
                    if(!(Meshfile >> Y_host(i)))
                    {
                        throw runtime_error("Y Coordinate reading failed");
                    }
                }
                Kokkos::deep_copy(X, X_host);
                Kokkos::deep_copy(Y, Y_host);

                imported.X = X;
                imported.Y = Y; 
                MeshData = imported;
            }
            else
            {
                throw runtime_error("Error: Mesh file could not be opened");
            }
            Meshfile.close();
        }
        catch (const exception& Error)
        {
            cerr << Error.what() << endl;
            throw;
        }
    }
    
    // Export physical cells only; preserve structured connectivity and scalar precision.
    template <typename T, typename CellData, typename InitData>
    void Export_Solution(const string& path, const MeshData_Structured<T>& meshdata, const CellData& cellsdata, const InitData& init_data, int iteration, double physical_time, double time_step)
    {
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "VTK export requires float or double");
        if (meshdata.Ni < 2 || meshdata.Nj < 2 || cellsdata.Ni_cells != meshdata.Ni - 1 || cellsdata.Nj_cells != meshdata.Nj - 1 || iteration < 0)
            throw std::invalid_argument("Solution export requires matching mesh/cell dimensions and a nonnegative iteration");
        const auto destination = std::filesystem::path(path);
        if (!destination.parent_path().empty())
            std::filesystem::create_directories(destination.parent_path());
        std::ofstream output(destination, std::ios::binary | std::ios::trunc);
        output.exceptions(std::ios::failbit | std::ios::badbit);
        output.imbue(std::locale::classic());
        output << std::setprecision(std::numeric_limits<T>::max_digits10);

        const auto x = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, meshdata.X);
        const auto y = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, meshdata.Y);
        const auto rho = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, cellsdata.rho_old);
        const auto momentum_u = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, cellsdata.Momentum_u_old);
        const auto momentum_v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, cellsdata.Momentum_v_old);
        const auto energy = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, cellsdata.E_old);
        const std::uint64_t nodes = std::uint64_t(meshdata.Ni) * std::uint64_t(meshdata.Nj);
        const std::uint64_t cells = std::uint64_t(cellsdata.Ni_cells) * std::uint64_t(cellsdata.Nj_cells);
        if (x.extent(0) != nodes || y.extent(0) != nodes)
            throw std::invalid_argument("Solution export found incomplete mesh coordinates");
        const char* vtk_type = sizeof(T) == 4 ? "Float32" : "Float64";
        const char* byte_order = std::endian::native == std::endian::little ? "LittleEndian" : "BigEndian";
        const T dynamic_pressure = T(0.5) * init_data.rho_inf * (init_data.velocity_u_inf * init_data.velocity_u_inf + init_data.velocity_v_inf * init_data.velocity_v_inf);
        const bool export_cp = dynamic_pressure > T{0} && std::isfinite(dynamic_pressure);
        const std::array<const char*, 9> names = {"density", "momentum_u", "momentum_v", "energy", "pressure", "sound_speed", "mach", "velocity", "cp"};
        const int array_count = export_cp ? 9 : 8;

        output << "<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"" << byte_order << "\" header_type=\"UInt64\">\n";
        output << "<StructuredGrid WholeExtent=\"0 " << meshdata.Ni - 1 << " 0 " << meshdata.Nj - 1 << " 0 0\">\n<FieldData>\n";
        const auto metadata = [&](const char* name, T value)
        {
            output << "<DataArray type=\"" << vtk_type << "\" Name=\"" << name << "\" NumberOfTuples=\"1\" format=\"ascii\">" << value << "</DataArray>\n";
        };
        metadata("gamma", init_data.gamma);
        metadata("mach_inf", init_data.Mach_inf);
        metadata("alpha_rad", init_data.alpha);
        metadata("rho_inf", init_data.rho_inf);
        metadata("pressure_inf", init_data.P_inf);
        if (physical_time >= 0 && std::isfinite(physical_time))
            output << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">" << std::setprecision(17) << physical_time << "</DataArray>\n";
        if (time_step >= 0 && std::isfinite(time_step))
            output << "<DataArray type=\"Float64\" Name=\"time_step\" NumberOfTuples=\"1\" format=\"ascii\">" << std::setprecision(17) << time_step << "</DataArray>\n";
        output << "<DataArray type=\"Int32\" Name=\"iteration\" NumberOfTuples=\"1\" format=\"ascii\">" << iteration << "</DataArray>\n</FieldData>\n";
        output << "<Piece Extent=\"0 " << meshdata.Ni - 1 << " 0 " << meshdata.Nj - 1 << " 0 0\">\n<PointData/>\n<CellData Scalars=\"pressure\" Vectors=\"velocity\">\n";
        std::uint64_t offset = 0;
        for (int field = 0; field < array_count; ++field)
        {
            const int components = field == 7 ? 3 : 1;
            output << "<DataArray type=\"" << vtk_type << "\" Name=\"" << names[field] << "\" NumberOfComponents=\"" << components << "\" format=\"appended\" offset=\"" << offset << "\"/>\n";
            offset += sizeof(std::uint64_t) + cells * components * sizeof(T);
        }
        output << "</CellData>\n<Points><DataArray type=\"" << vtk_type << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/></Points>\n</Piece>\n</StructuredGrid>\n<AppendedData encoding=\"raw\">_";

        // One reusable host buffer avoids an ASCII row per value and redundant connectivity.
        std::vector<T> buffer(std::size_t(3 * std::max(nodes, cells)));
        for (int field = 0; field < array_count; ++field)
        {
            const int components = field == 7 ? 3 : 1;
            for (int j = 0; j < cellsdata.Nj_cells; ++j)
                for (int i = 0; i < cellsdata.Ni_cells; ++i)
                {
                    const int k = i + cellsdata.N_ghostlayers + (cellsdata.Nitotal - 1) * (j + cellsdata.N_ghostlayers);
                    const std::size_t cell = std::size_t(i) + std::size_t(cellsdata.Ni_cells) * j;
                    const T density = rho(k), u = momentum_u(k) / density, v = momentum_v(k) / density;
                    const T pressure = (init_data.gamma - T{1}) * (energy(k) - T(0.5) * density * (u * u + v * v));
                    const T sound_speed = std::sqrt(init_data.gamma * pressure / density);
                    switch (field)
                    {
                    case 0: buffer[cell] = density; break;
                    case 1: buffer[cell] = momentum_u(k); break;
                    case 2: buffer[cell] = momentum_v(k); break;
                    case 3: buffer[cell] = energy(k); break;
                    case 4: buffer[cell] = pressure; break;
                    case 5: buffer[cell] = sound_speed; break;
                    case 6: buffer[cell] = std::sqrt(u * u + v * v) / sound_speed; break;
                    case 7: buffer[3 * cell] = u; buffer[3 * cell + 1] = v; buffer[3 * cell + 2] = T{0}; break;
                    case 8: buffer[cell] = (pressure - init_data.P_inf) / dynamic_pressure; break;
                    }
                }
            const std::uint64_t bytes = cells * components * sizeof(T);
            output.write(reinterpret_cast<const char*>(&bytes), sizeof(bytes));
            output.write(reinterpret_cast<const char*>(buffer.data()), std::streamsize(bytes));
        }
        for (std::uint64_t node = 0; node < nodes; ++node)
        {
            buffer[3 * node] = x(node);
            buffer[3 * node + 1] = y(node);
            buffer[3 * node + 2] = T{0};
        }
        const std::uint64_t point_bytes = nodes * 3 * sizeof(T);
        output.write(reinterpret_cast<const char*>(&point_bytes), sizeof(point_bytes));
        output.write(reinterpret_cast<const char*>(buffer.data()), std::streamsize(point_bytes));
        output << "</AppendedData>\n</VTKFile>\n";
        output.close();
    }


    template <typename T, typename InitData>
    SolutionSeriesWriter<T, InitData>::SolutionSeriesWriter(const std::filesystem::path& directory, const MeshData_Structured<T>& meshdata, const InitData& init_data) : directory(directory), meshdata(meshdata), init_data(init_data)
    {
        std::filesystem::create_directories(directory);
        if (std::filesystem::exists(directory / "solution.pvd"))
            throw std::runtime_error("A temporal solution already exists in the output directory; choose another directory");
        collection.exceptions(std::ios::failbit | std::ios::badbit);
        collection.open(directory / "solution.pvd", std::ios::in | std::ios::out | std::ios::trunc | std::ios::binary);
        collection.imbue(std::locale::classic());
        collection << std::setprecision(17) << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n";
        collection_end = collection.tellp();
        collection << "</Collection>\n</VTKFile>\n" << std::flush;
    }

    template <typename T, typename InitData>
    template <typename CellData>
    void SolutionSeriesWriter<T, InitData>::write(const CellData& cellsdata, int iteration, double physical_time, double time_step)
    {
        if (iteration <= last_iteration || !std::isfinite(physical_time) || !std::isfinite(time_step) || physical_time < 0 || time_step < 0 || (last_iteration >= 0 && physical_time <= last_time))
            throw std::invalid_argument("Temporal samples require increasing iterations and physical times, and finite nonnegative time steps");
        std::ostringstream filename;
        filename << "solution_" << std::setfill('0') << std::setw(8) << iteration << ".vts";
        const auto destination = directory / filename.str();
        if (std::filesystem::exists(destination))
            throw std::runtime_error("A temporal sample already exists; choose another output directory");
        Export_Solution<T>(destination.string(), meshdata, cellsdata, init_data, iteration, physical_time, time_step);
        collection.seekp(collection_end);
        collection << "<DataSet timestep=\"" << physical_time << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>\n";
        collection_end = collection.tellp();
        collection << "</Collection>\n</VTKFile>\n" << std::flush;
        last_time = physical_time;
        last_iteration = iteration;
    }

}