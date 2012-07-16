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


#ifndef TESTSPHEROIDEXPERIMENTS_HPP_
#define TESTSPHEROIDEXPERIMENTS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "OffLatticeSimulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "AveragedSourcePde.hpp"
#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellBasedPdeHandler.hpp"
#include "ChastePoint.hpp"
#include "SmartPointers.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"

class TestSpheroidExperiments : public AbstractCellBasedTestSuite
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

    void TestMeshBasedSpheroidWithPde() throw(Exception)
    {
        // Create a simple 3D mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);

        // Set up cells
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            StochasticOxygenBasedCellCycleModel* p_model = new StochasticOxygenBasedCellCycleModel();
            p_model->SetDimension(3);
            p_model->SetStemCellG1Duration(2.0);
            p_model->SetHypoxicConcentration(0.1);
            p_model->SetQuiescentConcentration(0.3);
            p_model->SetCriticalHypoxicDuration(8);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (  p_model->GetStemCellG1Duration()
                                 + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);
        cell_population.SetOutputCellVolumes(true);
        cell_population.SetWriteVtkAsPoints(true);
        
        // Set up cell data on the cell population
        cell_population.SetDataOnAllCells("oxygen", 1.0);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetEndTime(100);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetOutputDirectory("Plos2012_MeshBasedSpheroidWithPde");

        // Set up PDE and pass to simulation via handler
        CellwiseSourcePde<3> pde(cell_population, -1.0);
        ConstBoundaryCondition<3> bc(1.0);
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<3> pde_and_bc(&pde, &bc, is_neumann_bc);
        pde_and_bc.SetDependentVariableName("oxygen");

        CellBasedPdeHandler<3> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<3>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Solve the system
        simulator.Solve();

        // Save results
        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);

        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();
    }
};

#endif /*TESTSPHEROIDEXPERIMENTS_HPP_*/
