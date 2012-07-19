/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#ifndef TESTELECTROMECHANICSLITERATEPAPER_HPP_
#define TESTELECTROMECHANICSLITERATEPAPER_HPP_


#include <cxxtest/TestSuite.h>
#include "PlaneStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "FileFinder.hpp"
#include "VtkMeshWriter.hpp"

/*
* This example is based on one of the cardiac electro-mechanics tutorials,
* see these for further details.
*
* Remember to run with build=GccOptNative_ndebug for speed
*/
class TestElectroMechanicsLiteratePaper : public CxxTest::TestSuite
{
public:
    void TestTwistingCube() throw(Exception)
    {
    	// This class sets up cells with the correct stimuli to
        // stimulate just those cells on the x=0 surface at t=0.
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3u> cell_factory(-1000*1000);

        // Set up two meshes of 1mm by 1mm by 1mm
        // We have a fine electrics mesh and a coarser mechanics one.
        TetrahedralMesh<3u,3u> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01, 0.1, 0.1, 0.1);

        QuadraticMesh<3> mechanics_mesh(0.02, 0.1, 0.1, 0.1);

        // fix the nodes on Z=0
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<3u>::GetNodesByComponentValue(mechanics_mesh, 2, 0.0);

        HeartConfig::Instance()->SetSimulationDuration(36.0); //ms

        // Define the electro-mechanics problem
        ElectroMechanicsProblemDefinition<3u> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0); // ms

        std::string output_directory = "Plos2012_ElectroMechanics";
        std::string fibre_file_name = "5by5by5_fibres.ortho";

        // This is how to generate a fibre file for this mesh
        // (this could be generated once, stored and re-loaded as below)
        // We use a Streeter-style formula here.
        {
            OutputFileHandler handler(output_directory + "Fibres");
            out_stream p_file = handler.OpenOutputFile(fibre_file_name);
            *p_file << mechanics_mesh.GetNumElements() << "\n"; // first line is number of entries
            std::vector<c_vector<double,3u> > fibre_directions;
            for(unsigned i=0; i<mechanics_mesh.GetNumElements(); i++)
            {
                double X = mechanics_mesh.GetElement(i)->CalculateCentroid()(0);
                double theta = M_PI/3 - 10*X*2*M_PI/3; // 60 degrees when X=0, -60 when X=0.1;

                c_vector<double,3u> fibre_direction;
                fibre_direction[0] = 0;
                fibre_direction[1] = cos(theta);
                fibre_direction[2] = sin(theta);
                *p_file <<  fibre_direction[0] << " " << fibre_direction[1]  << " " << fibre_direction[2]  // first three entries are fibre direction
                        << " 0 " << -fibre_direction[2] << " " << fibre_direction[1]                       // next three are sheet direction
                        << " 1 0 0\n";                                                                     // then normal to sheet direction
                fibre_directions.push_back(fibre_direction);
            }
            p_file->close();

            // We only compile the following if VTK is installed and set up.
            // This is optional - just for visualizing the fibre directions as in the paper.
            // the simulation will run without it.
#ifdef CHASTE_VTK
            VtkMeshWriter<3u,3u> mesh_writer(output_directory+ "Fibres", "mesh", false);

            mesh_writer.AddCellData("Fibre Directions", fibre_directions);
            mesh_writer.WriteFilesUsingMesh(mechanics_mesh);
#endif // CHASTE_VTK
        }

        // Load up the file we just wrote to use as fibre directions for this mechanics problem.
        FileFinder fibre_file_finder(output_directory + "Fibres/" + fibre_file_name, RelativeTo::ChasteTestOutput);
        problem_defn.SetVariableFibreSheetDirectionsFile(fibre_file_finder, false);

        // Set up and solve the full cardiac electro-mechanics problem.
        CardiacElectroMechanicsProblem<3u> problem(COMPRESSIBLE,
                                                   MONODOMAIN,
                                                   &electrics_mesh,
                                                   &mechanics_mesh,
                                                   &cell_factory,
                                                   &problem_defn,
                                                   output_directory);

        // Run the simulation
        problem.Solve();

        // Report where time was spent to std::cout.
        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }
};

#endif /* TESTELECTROMECHANICSLITERATEPAPER_HPP_ */
