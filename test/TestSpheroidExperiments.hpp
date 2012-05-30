#ifndef TESTSPHEROID_HPP_
#define TESTSPHEROID_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"


#include "OffLatticeSimulation.hpp"

#include "CellBasedEventHandler.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellwiseSourcePde.hpp"
#include "AveragedSourcePde.hpp"
#include "AveragedPointSourcePde.hpp"
#include "CellData.hpp"
#include "VolumeDependentAveragedSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "WildTypeCellMutationState.hpp"
#include "MutableMesh.hpp"
#include "LogFile.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedPdeHandler.hpp"
#include "ChasteCuboid.hpp"
#include "ChastePoint.hpp"
#include "SmartPointers.hpp"
#include "OxygenBasedCellKiller.hpp"

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
		// Create mesh
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> start_mesh;
        start_mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("TestSolveMethodSpheroidSimulation3DMesh", "StartMesh");
        mesh_writer.WriteFilesUsingMesh(start_mesh);

		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(start_mesh);

		// Set up cells
		std::vector<CellPtr> cells;
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
			p_model->SetDimension(3);
			p_model->SetCellProliferativeType(STEM);
			p_model->SetStemCellG1Duration(2.0);

			CellPtr p_cell(new Cell(p_state, p_model));
			double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}

		// Set up cell population
		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		cell_population.SetMechanicsCutOffLength(1.5);
		cell_population.SetAbsoluteMovementThreshold(DBL_MAX);
		// Set up cell data on the cell population
		MAKE_PTR_ARGS(CellData, p_cell_data, (1));
		p_cell_data->SetItem(0, 1.0);
		cell_population.AddClonedDataToAllCells(p_cell_data);

		// Set up PDE
		AveragedSourcePde<3> pde(cell_population, -1.0);
		ConstBoundaryCondition<3> bc(1.0);
		bool is_neumann_bc = false;
		PdeAndBoundaryConditions<3> pde_and_bc(&pde, &bc, is_neumann_bc);

		std::vector<PdeAndBoundaryConditions<3>*> pde_and_bc_collection;
		pde_and_bc_collection.push_back(&pde_and_bc);

        CellBasedPdeHandler<3> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);

		// Set up cell-based simulation
		OffLatticeSimulation<3> simulator(cell_population);
		simulator.SetEndTime(100.0);
		simulator.SetSamplingTimestepMultiple(120);
		simulator.SetCellBasedPdeHandler(&pde_handler);

		// Set output directory
		simulator.SetOutputDirectory("Simple3DMonolayerExperiment");


        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

		// Tell simulator to use the coarse mesh.
        ChastePoint<3> lower(0.0, 0.0, 0.0);
        ChastePoint<3> upper(50.0, 50.0, 50.0);
        ChasteCuboid<3> cuboid(lower, upper);

		pde_handler.UseCoarsePdeMesh(10.0, cuboid, true);

		//Solve the system
		simulator.Solve();

		CellBasedEventHandler::Headings();
		CellBasedEventHandler::Report();
    }
};

#endif /*TESTSPHEROID_HPP_*/
