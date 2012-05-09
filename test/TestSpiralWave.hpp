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


#ifndef TESTSPIRALWAVE_HPP_
#define TESTSPIRALWAVE_HPP_

#include <cxxtest/TestSuite.h>

#include "MonodomainProblem.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "LuoRudyCellFactory.hpp" // This is defined in this project 'src' folder, the rest are in Chaste 3.0.
#include "PetscSetupAndFinalize.hpp"

/*
 * Having included all the necessary header files, we proceed by defining the test class.
 */
class TestSpiralWave : public CxxTest::TestSuite
{
public:
    void TestSpiralWaveSimulation() throw (Exception)
    {
        /*
         * We will auto-generate a mesh this time, and pass it in, rather than
         * provide a mesh file name. This is how to generate a cuboid mesh with
         * a given spatial stepsize h
         */
        DistributedTetrahedralMesh<2,2> mesh;
        double node_spacing_in_mesh = 0.015;
        double mesh_width = 3; // cm
        mesh.ConstructRegularSlabMesh(node_spacing_in_mesh, mesh_width /*length*/, mesh_width /*width*/);
        /*
         * Set the simulation duration, etc, and create an instance of the cell factory.
         * One thing that should be noted for monodomain problems, the ''intracellular
         * conductivity'' is used as the monodomain effective conductivity (not a
         * harmonic mean of intra and extracellular conductivities).
         * So if you want to alter the monodomain conductivity call
         * `HeartConfig::Instance()->SetIntracellularConductivities`
         */
        HeartConfig::Instance()->SetSimulationDuration(500); //ms
        HeartConfig::Instance()->SetOutputDirectory("Plos2012_SpiralWave");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 1);

        LuoRudyCellFactory cell_factory(mesh_width,mesh_width);

        /*
         * Now we declare the problem class, `MonodomainProblem<2>` instead of `BidomainProblem<2>`.
         * The interface for both is the same.
         */
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        /*
         * If a mesh-file-name hasn't been set using `HeartConfig`, we have to pass in
         * a mesh using the `SetMesh` method (must be called before `Initialise`).
         */
        monodomain_problem.SetMesh(&mesh);

//        /*
//         * By default data for all nodes is output, but for big simulations, sometimes this
//         * might not be required, and the action potential only at certain nodes required.
//         * The following code shows how to output the results at the first, middle and last
//         * nodes, for example. (The output is written to the HDF5 file; no meshalyzer output
//         * will be made. HDF5 files can be read using Matlab). We are not using this in this
//         * simulation however (hence the boolean being set to false).
//         */
//        bool partial_output = false;
//        if(partial_output)
//        {
//            std::vector<unsigned> nodes_to_be_output;
//            nodes_to_be_output.push_back(0);
//            nodes_to_be_output.push_back((unsigned)round( (mesh.GetNumNodes()-1)/2 ));
//            nodes_to_be_output.push_back(mesh.GetNumNodes()-1);
//            monodomain_problem.SetOutputNodes(nodes_to_be_output);
//        }

        /* `SetWriteInfo` is a useful method that means that the min/max voltage is
         * printed as the simulation runs (useful for verifying that cells are stimulated
         * and the wave propagating, for example) (although note scons does buffer output
         * before printing to screen) */
        monodomain_problem.SetWriteInfo();

        /* Finally, call `Initialise` and `Solve` as before */
        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        /* This part is just to check nothing has accidentally been changed in this example */
        ReplicatableVector voltage(monodomain_problem.GetSolution());
        TS_ASSERT_DELTA(voltage[0], 35.0939, 1e-2);
    }
};

#endif /*TESTSPIRALWAVE_HPP_*/
