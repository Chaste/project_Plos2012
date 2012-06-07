#ifndef TESTSPHEROIDEXPERIMENTS_HPP_
#define TESTSPHEROIDEXPERIMENTS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "OffLatticeSimulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "AveragedSourcePde.hpp"
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

	void TestSpheroidWithPde() throw(Exception)
    {
        // Create a single node
        std::vector<Node<3>* > nodes;
        nodes.push_back(new Node<3>(0, false,  0.0, 0.0, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

		// Create a corresponding cell
		std::vector<CellPtr> cells;

		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		StochasticOxygenBasedCellCycleModel* p_model = new StochasticOxygenBasedCellCycleModel();
		p_model->SetDimension(3);
		p_model->SetCellProliferativeType(STEM);
		p_model->SetStemCellG1Duration(2.0);

		CellPtr p_cell(new Cell(p_state, p_model));
		double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
		p_cell->SetBirthTime(birth_time);
		cells.push_back(p_cell);

		// Set up cell population
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.SetMechanicsCutOffLength(1.5);
		cell_population.SetAbsoluteMovementThreshold(DBL_MAX);

        // Set up cell data on the cell population
        cell_population.SetDataOnAllCells("oxygen", 1.0);

		// Set up cell-based simulation
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.SetEndTime(100);
		simulator.SetSamplingTimestepMultiple(120);

        // Set up PDE and pass to simulation via handler
        AveragedSourcePde<3> pde(cell_population, -1.0);
        ConstBoundaryCondition<3> bc(1.0);
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<3> pde_and_bc(&pde, &bc, is_neumann_bc);
        pde_and_bc.SetDependentVariableName("oxygen");

        CellBasedPdeHandler<3> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
		simulator.SetCellBasedPdeHandler(&pde_handler);

		// Set output directory
		simulator.SetOutputDirectory("TestSpheroidWithPde");

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

		// Use a coarse mesh when solving the PDE
        ChastePoint<3> lower(0.0, 0.0, 0.0);
        ChastePoint<3> upper(50.0, 50.0, 50.0);
        ChasteCuboid<3> cuboid(lower, upper);
		pde_handler.UseCoarsePdeMesh(10.0, cuboid, true);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<3>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

		// Solve the system
		simulator.Solve();

		CellBasedEventHandler::Headings();
		CellBasedEventHandler::Report();
    }
};

#endif /*TESTSPHEROIDEXPERIMENTS_HPP_*/
