/*
 * VTKWriter.cpp
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#include "VTKWriter.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <spdlog/spdlog.h>

namespace outputWriter {

VTKWriter::VTKWriter() = default;

VTKWriter::~VTKWriter() = default;

VTKFile_t VTKWriter::initializeOutput(int numParticles) {

    VTKFile_t vtk_file("UnstructuredGrid");

    // per point, we add type, position, velocity and force
    PointData pointData;
    DataArray_t mass(type::Float32, "mass", 1);
    DataArray_t velocity(type::Float32, "velocity", 3);
    DataArray_t forces(type::Float32, "force", 3);
    DataArray_t type(type::Int32, "type", 1);
    pointData.DataArray().push_back(mass);
    pointData.DataArray().push_back(velocity);
    pointData.DataArray().push_back(forces);
    pointData.DataArray().push_back(type);

    CellData cellData; // we don't have cell data => leave it empty

    // 3 coordinates
    Points points;
    DataArray_t pointCoordinates(type::Float32, "points", 3);
    points.DataArray().push_back(pointCoordinates);

    Cells cells; // we don't have cells, => leave it empty
    // for some reasons, we have to add a dummy entry for paraview
    DataArray_t cells_data(type::Float32, "types", 0);
    cells.DataArray().push_back(cells_data);

    PieceUnstructuredGrid_t piece(pointData, cellData, points, cells,
                                  numParticles, 0);
    UnstructuredGrid_t unstructuredGrid(piece);
    vtk_file.UnstructuredGrid(unstructuredGrid);

    return vtk_file;
  }

  void VTKWriter::writeFile(const std::string &filename, int iteration, std::vector<Particle> &particles) {
      std::stringstream strstr;
      strstr << filename << "_" << std::setfill('0') << std::setw(4) << iteration << ".vtu";

      int numParticles = 0;
      for (const auto& particle : particles) {
          if (particle.isInDomain()) numParticles++;
      }
      auto vtkFile = initializeOutput(numParticles);

      for(Particle& p: particles) {
          if (p.isInDomain()) {
              plotParticle(vtkFile, p);
          }
      }

      std::ofstream file(strstr.str().c_str());
      VTKFile(file, vtkFile);
      file.close();

    //delete &vtkFile;

  }

  void VTKWriter::plotParticle(VTKFile_t& vtk_file, Particle &p) {
    if (vtk_file.UnstructuredGrid().present()) {
      spdlog::trace("UnstructuredGrid is present");
    } else {
      spdlog::error("ERROR: No UnstructuredGrid present");
    }

    PointData::DataArray_sequence &pointDataSequence =
        vtk_file.UnstructuredGrid()->Piece().PointData().DataArray();
    PointData::DataArray_iterator dataIterator = pointDataSequence.begin();

    dataIterator->push_back(p.getM());
    // cout << "Appended mass data in: " << dataIterator->Name();

    dataIterator++;
    dataIterator->push_back(p.getV()[0]);
    dataIterator->push_back(p.getV()[1]);
    dataIterator->push_back(p.getV()[2]);
    // cout << "Appended velocity data in: " << dataIterator->Name();

    dataIterator++;
    dataIterator->push_back(p.getOldF()[0]);
    dataIterator->push_back(p.getOldF()[1]);
    dataIterator->push_back(p.getOldF()[2]);
    // cout << "Appended force data in: " << dataIterator->Name();

    dataIterator++;
    dataIterator->push_back(p.getType());

    Points::DataArray_sequence &pointsSequence =
        vtk_file.UnstructuredGrid()->Piece().Points().DataArray();
    Points::DataArray_iterator pointsIterator = pointsSequence.begin();
    pointsIterator->push_back(p.getX()[0]);
    pointsIterator->push_back(p.getX()[1]);
    pointsIterator->push_back(p.getX()[2]);
    // if(p.getV()[1] < -100000){
    // spdlog::warn("particle position: {} {} {}, velocity: {} {} {}, force: {} {} {}",
    //         (p.getX()[0]),(p.getX()[1]),(p.getX()[2]),
    //         (p.getV()[0]),(p.getV()[1]),(p.getV()[2]),
    //         (p.getF()[0]),(p.getF()[1]),(p.getF()[2]));
    // }

  }

  } // namespace outputWriter
