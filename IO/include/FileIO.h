#pragma once

#include <iostream>
#include <filesystem>
#include <fstream>
#include <MeshHandler.h>

using namespace std;
using namespace Mesh;
namespace IO
{
    std::filesystem::path Create_RunDirectory(const std::filesystem::path& root);

    template <typename T, MeshFileFormat Format = MeshFileFormat::Plot3D>
    void Import_Mesh(const string& path, MeshData_Structured<T>& MeshData);

    template <typename T, typename CellData, typename InitData>
    void Export_Solution(const string& path, const MeshData_Structured<T>& meshdata, const CellData& cellsdata, const InitData& init_data, int iteration, double physical_time = -1, double time_step = -1);

    template <typename T, typename InitData>
    class SolutionSeriesWriter
    {
        public:
            SolutionSeriesWriter(const std::filesystem::path& directory, const MeshData_Structured<T>& meshdata, const InitData& init_data);
            template <typename CellData>
            void write(const CellData& cellsdata, int iteration, double physical_time, double time_step);
        private:
            std::filesystem::path directory;
            const MeshData_Structured<T>& meshdata;
            const InitData& init_data;
            std::fstream collection;
            std::streampos collection_end;
            double last_time = -1;
            int last_iteration = -1;
    };

}
#include "FileIO.tpp"
