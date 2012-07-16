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


#ifndef TESTCELLBASEDMULTIPLECRYPT_HPP_
#define TESTCELLBASEDMULTIPLECRYPT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "SimplifiedDeltaNotchOffLatticeSimulation.hpp"
#include "DeltaNotchOffLatticeSimulation.hpp"
#include "CryptCellsGenerator.hpp"
#include "CellsGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "SimpleWntCellCycleModelWithDeltaNotch.hpp"
#include "DeltaNotchCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "MultipleCryptGeometryBoundaryCondition.hpp"

#include "Debug.hpp"

class TestCellBasedMultipleCrypt : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void Test3dCrypt() throw (Exception)
    {
        double crypt_length = 4.0;
        double crypt_radius = 1.0;
        double villus_length = 10.0;
        double villus_radius = 2.0;
        double domain_width = 12.0;
        double domain_height = 2.0*crypt_radius+crypt_length+villus_length+2.0*villus_radius;

        // Put a single cell at the base of each crypt. TODO Change this back to one cell and check its OK
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u,  false,  0.5*domain_width, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  0.5*domain_width+0.1, 0.1, 0.0));
        nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5*domain_width, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  0.1, 0.5*domain_width+0.1, 0.0));
        nodes.push_back(new Node<3>(4u,  false,  0.5*domain_width, domain_width, 0.0));
        nodes.push_back(new Node<3>(5u,  false,  0.5*domain_width+0.1, domain_width-0.1, 0.0));
        nodes.push_back(new Node<3>(6u,  false,  domain_width, 0.5*domain_width, 0.0));
        nodes.push_back(new Node<3>(7u,  false,  domain_width-0.1, 0.5*domain_width+0.1, 0.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            SimpleWntCellCycleModelWithDeltaNotch* p_model = new SimpleWntCellCycleModelWithDeltaNotch();
            p_model->SetDimension(3);

            /* We choose to initialise the concentrations to random levels in each cell 
               TODO check this does work with random ICS. */
            std::vector<double> initial_conditions;
            initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            p_model->SetInitialConditions(initial_conditions);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            p_cell->SetCellProliferativeType(TRANSIT);
            cells.push_back(p_cell);
        }

        // Create a node-based cell population
        NodeBasedCellPopulation<3> crypt(mesh, cells);
        crypt.SetMechanicsCutOffLength(1.5);

        crypt.SetAbsoluteMovementThreshold(10);

        // Output some useful information for plotting in VTK format.
        crypt.SetOutputCellProliferativeTypes(true);
        crypt.SetOutputCellMutationStates(true);
        crypt.SetOutputCellAncestors(true);


        // Set up cell-based simulation
		SimplifiedDeltaNotchOffLatticeSimulation<3> simulator(crypt);
        simulator.SetOutputDirectory("Plos2012_MultipleCrypt");
        simulator.SetDt(1.0/100.0);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetOutputNodeVelocities(true);

        // Create a force law and pass it to the simulation
        // We use linear springs between cells up to 1.5 (cell diameters) apart
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(30.0); //normally 15.0;
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        MAKE_PTR_ARGS(MultipleCryptGeometryBoundaryCondition, 
                      p_boundary_condition, 
                      (&crypt, crypt_radius, crypt_length, villus_radius, villus_length, domain_width));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);


        // Main sink of cells is at the top of the villus
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, 
                      p_cell_killer_1,
                      (&crypt, (domain_height-0.25)*unit_vector<double>(3,2), unit_vector<double>(3,2)));
        simulator.AddCellKiller(p_cell_killer_1);

        // Create an instance of a Wnt concentration, this dictates where cell division occurs
        // in the SimpleWntCellCycleModelWithDeltaNotch
        WntConcentration<3>::Instance()->SetType(LINEAR);
        WntConcentration<3>::Instance()->SetCellPopulation(crypt);
        WntConcentration<3>::Instance()->SetCryptLength(crypt_length+2.0*crypt_radius);

        // Run simulation
        simulator.SetEndTime(250.0);
        simulator.Solve(); // to 250 hours

        // Add a random cell killer to represent random death in the epithelial layer.
        MAKE_PTR_ARGS(RandomCellKiller<3>, p_cell_killer_2,(&crypt, 0.005)); // prob of death in an hour
        simulator.AddCellKiller(p_cell_killer_2);

        // Label each cell according to its current node index so we can track clonal spread
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        simulator.SetEndTime(1000.0);
        simulator.Solve(); // to 1000 hours
    }

};

#endif /*TESTCELLBASEDMULTIPLECRYPT_HPP_*/
