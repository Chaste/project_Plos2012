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

#ifndef TESTDELTANOTCHSIMULATION_HPP_
#define TESTDELTANOTCHSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SmartPointers.hpp"

/*
 * The next header file defines a simple stochastic cell-cycle model that includes the functionality
 * for solving each cell's Delta/Notch signalling ODE system at each time step, using information about neighbouring
 * cells through the {{{CellwiseData}}} singleton. We note that in this simple cell-cycle model, the
 * proliferative status of each cell is unaffected by its Delta/Notch activity; such dependence could
 * easily be introduced given an appropriate model of this coupling.
 */
#include "DeltaNotchCellCycleModel.hpp"
/*
 * The next header file defines the class that simulates the evolution of a {{{CellPopulation}}},
 * specialized to deal with updating of the {{{CellwiseData}}} singleton to deal with Delta-Notch
 * signalling between cells.
 */
#include "DeltaNotchOffLatticeSimulation.hpp"

/*
 * Having included all the necessary header files, we proceed by defining the test class.
 */
class TestDeltaNotchSimulation : public AbstractCellBasedTestSuite
{
public:
    /*
     * == A node-based monolayer with Delta/Notch signalling ==
     *
     * In the next test we run a similar simulation as the delta-notch tutorial,
     * here with node-based 'overlapping spheres' model.
     */
    void TestNodeBasedMonolayerWithDeltaNotch() throw (Exception)
    {
        /*
         * Most of the code in this test is the same as in the previous test,
         * except we now create a 'nodes-only mesh' and {{{NodeBasedCellPopulation}}}.
         */
        HoneycombMeshGenerator generator(5, 5);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        /* The mechanics cut-off length is also used in this simulation to determine nearest
         * neighbours for the purpose of the Delta/Notch intercellular signalling model.
         */
        cell_population.SetMechanicsCutOffLength(1.5);

        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAges(true);

        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(p_mesh->GetNumNodes(), 3);
        p_data->SetCellPopulation(&cell_population);
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 0);
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 1);
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 2);
        }

        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("DeltaNotchSimulation");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(50.0);

        /* As we are using a node-based cell population, we use an appropriate force law. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        simulator.Solve();

        /* Finally, as before, we call {{{Destroy()}}} on any singleton classes. */
        CellwiseData<2>::Destroy();

        /* To avoid memory leaks, we also delete any pointers we created in the test. */
        delete p_mesh;
    }
    /*
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information
     *
     * Load the file {{{/tmp/$USER/testoutput/TestNodeBasedMonolayerWithDeltaNotch/results_from_time_0/results.pvd}}},
     * add a spherical glyph.
     */
};

#endif /*TESTDELTANOTCHSIMULATION_HPP_*/
