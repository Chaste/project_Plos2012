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

#include <iomanip>
#include <boost/foreach.hpp>

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

    /*
     * These methods are cxx-test instructions
     * running before and after each test below.
     * They just report the time the
     * test took.
     */
    void setUp()
    {
        AbstractCellBasedTestSuite::setUp();
        CellBasedEventHandler::Reset();
    }
    void tearDown()
    {
        AbstractCellBasedTestSuite::tearDown();

        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();
    }

public:

    /*
     * This example is split into two separate codes to showcase the checkpointing
     * abilities of Chaste.
     *
     * The first test runs from t=0 to t=100,
     * and the second from t=100 to t=160.
     *
     * It could equally well be reproduced by setting the end time in the first
     * test to 160.
     */
    void TestMeshBasedSpheroidWithPde() throw(Exception)
    {
        // Create a simple 3D mesh, initially comprised of just five nodes
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);

        // Set up cells for each of the nodes in the mesh
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

        // Create cell population - a mapping between a mesh and cells.
        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);
        cell_population.SetOutputCellVolumes(true);
        cell_population.SetWriteVtkAsPoints(true);
        
        // Set up cell data on the cell population
        cell_population.SetDataOnAllCells("oxygen", 1.0);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetEndTime(100); // hours
        simulator.SetSamplingTimestepMultiple(120); // Default timestep is 30 seconds,
                                                    //so this gives one set of output each hour.
        simulator.SetOutputDirectory("Plos2012_MeshBasedSpheroidWithPde");

        // Set up PDE and boundary conditions
        CellwiseSourcePde<3> pde(cell_population, -1.0);
        ConstBoundaryCondition<3> bc(1.0);
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<3> pde_and_bc(&pde, &bc, is_neumann_bc);
        pde_and_bc.SetDependentVariableName("oxygen");

        // Create a handler (for any number of PDEs+BCs, in this case we just add one).
        CellBasedPdeHandler<3> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Pass PDE handler to the simulation
        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<3>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run the simulation
        simulator.Solve();

        // Save results
        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);
    }

    void TestLongerMeshBasedSpheroidWithPde() throw(Exception)
    {
        // The archive is to be copied from previous test output
        // It could be stored and re-loaded from anywhere you like
        // if you want to experiment with different interventions on
        // an existing spheroid.
        FileFinder test_data_directory("Plos2012_MeshBasedSpheroidWithPde/archive",
                                       RelativeTo::ChasteTestOutput);

        // to the testoutput/archive directory to continue running the simulation
        OutputFileHandler archive_handler("Plos2012_LongerMeshBasedSpheroidWithPde/archive");

        // Following is done in two lines to avoid a bug in Intel compiler v12.0!
        std::vector<FileFinder> temp_files = test_data_directory.FindMatches("*");
        BOOST_FOREACH(FileFinder temp_file, temp_files)
        {
            archive_handler.CopyFileTo(temp_file);
        }

        // Load the simulation up
        OffLatticeSimulation<3>* p_simulator
            = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load("Plos2012_LongerMeshBasedSpheroidWithPde", 100);

        // Change some settings
        p_simulator->SetEndTime(160);
        p_simulator->SetOutputDirectory("Plos2012_LongerMeshBasedSpheroidWithPde");

        // Solve the system
        p_simulator->Solve();

        // Save results
        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulator);
    }
};

#endif /*TESTSPHEROIDEXPERIMENTS_HPP_*/
