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
        double villus_length = 8.0;
        double villus_radius = 2.0;
        double domain_width = 12.0;
        double domain_height = 2.0*crypt_radius+crypt_length+villus_length+2.0*villus_radius;

        // Put a single cell at the base of each crypt.
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u,  false,  0.5*domain_width, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  0.0, 0.5*domain_width, 0.0));
        nodes.push_back(new Node<3>(2u,  false,  0.5*domain_width, domain_width, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  domain_width, 0.5*domain_width, 0.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            SimpleWntCellCycleModelWithDeltaNotch* p_model = new SimpleWntCellCycleModelWithDeltaNotch();
            p_model->SetCellProliferativeType(TRANSIT);
            p_model->SetDimension(3);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        MAKE_PTR(CellLabel, p_label);

        // Create a node-based cell population
        NodeBasedCellPopulation<3> crypt(mesh, cells);
        crypt.SetMechanicsCutOffLength(1.5);
        crypt.SetOutputCellProliferativeTypes(true);
        crypt.SetOutputCellMutationStates(true);
        crypt.SetOutputCellAncestors(true);
        crypt.SetAbsoluteMovementThreshold(10);

        /* We choose to initialise the concentrations to random levels in each cell. */
		 for (AbstractCellPopulation<3>::Iterator cell_iter = crypt.Begin();
			  cell_iter != crypt.End();
			  ++cell_iter)
		 {
			 cell_iter->GetCellData()->SetItem("notch", RandomNumberGenerator::Instance()->ranf());
			 cell_iter->GetCellData()->SetItem("delta", RandomNumberGenerator::Instance()->ranf());
			 cell_iter->GetCellData()->SetItem("mean delta", RandomNumberGenerator::Instance()->ranf());
		 }


        // Set up cell-based simulation
		SimplifiedDeltaNotchOffLatticeSimulation<3> simulator(crypt);
        simulator.SetOutputDirectory("Plos2012_MultipleCrypt");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetOutputNodeVelocities(true);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(30.0); //normally 15.0;
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        MAKE_PTR_ARGS(MultipleCryptGeometryBoundaryCondition, p_boundary_condition, (&crypt, crypt_radius, crypt_length, villus_radius, villus_length, domain_width));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);


        // Main sink of cells is at the top of the villus
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer_1,(&crypt, (domain_height-0.25)*unit_vector<double>(3,2), unit_vector<double>(3,2)));
        simulator.AddCellKiller(p_cell_killer_1);



        // Create an instance of a Wnt concentration
        WntConcentration<3>::Instance()->SetType(LINEAR);
        WntConcentration<3>::Instance()->SetCellPopulation(crypt);
        WntConcentration<3>::Instance()->SetCryptLength(crypt_length+2*crypt_radius);

        // Run simulation
        simulator.SetEndTime(100);
        simulator.Solve(); // to 200

        // Label each cell according to its current node index
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();
        simulator.SetEndTime(1000.0);

        MAKE_PTR_ARGS(RandomCellKiller<3>, p_cell_killer_2,(&crypt, 0.008)); // prob of death in an hour
        simulator.AddCellKiller(p_cell_killer_2);

        // Run simulation
        simulator.Solve();
    }

};

#endif /*TESTCELLBASEDMULTIPLECRYPT_HPP_*/
