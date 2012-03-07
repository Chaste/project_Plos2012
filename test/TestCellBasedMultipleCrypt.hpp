#ifndef TESTCELLBASEDMULTIPLECRYPT_HPP_
#define TESTCELLBASEDMULTIPLECRYPT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "CellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "SimpleWntCellCycleModel.hpp"
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
        double crypt_length = 1.0;
        double crypt_radius = 1.0;

        // Create mesh Need so many nodes to ensure there are some node pairs setup
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u,  false,  2.0*crypt_radius, 2.0*crypt_radius, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  2.0*crypt_radius+1, 2.0*crypt_radius, 0.0));

        nodes.push_back(new Node<3>(2u,  false,  2.0*crypt_radius, 10.0*crypt_radius, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  2.0*crypt_radius+1.0, 10.0*crypt_radius, 0.0));


        nodes.push_back(new Node<3>(4u,  false,  10.0*crypt_radius, 2.0*crypt_radius, 0.0));
        nodes.push_back(new Node<3>(5u,  false,  10.0*crypt_radius+1.0, 2.0*crypt_radius, 0.0));

        nodes.push_back(new Node<3>(6u,  false,  10.0*crypt_radius, 10.0*crypt_radius, 0.0));
        nodes.push_back(new Node<3>(7u,  false,  10.0*crypt_radius+1.0, 10.0*crypt_radius, 0.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<SimpleWntCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), TRANSIT);

        MAKE_PTR(CellLabel, p_label);

        // Create a node-based cell population
        NodeBasedCellPopulation<3> crypt(mesh, cells);
        crypt.SetMechanicsCutOffLength(1.5);
        crypt.SetOutputCellProliferativeTypes(true);
        crypt.SetOutputCellMutationStates(true);
        crypt.SetOutputCellAncestors(true);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(crypt);
        simulator.SetOutputDirectory("MultipleCryptDemo3d");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(200);
        simulator.SetOutputNodeVelocities(true);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(30.0); //normally 15.0;
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);


        MAKE_PTR_ARGS(MultipleCryptGeometryBoundaryCondition<3>, p_boundary_condition, (&crypt, crypt_radius, crypt_length));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        //        MAKE_PTR_ARGS(MultipleCryptGeometryBoundaryCondition<3>, p_boundary_condition, (&crypt, crypt_radius, crypt_length));
        //        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);
        //
        //        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition_1,(&crypt, 12*crypt_radius*unit_vector<double>(3,0), unit_vector<double>(3,0)));
        //        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_1);
        //
        //        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition_2,(&crypt, 12*crypt_radius*unit_vector<double>(3,1), unit_vector<double>(3,1)));
        //        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_2);
        //
        //        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition_3,(&crypt, zero_vector<double>(3), -unit_vector<double>(3,0)));
        //        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_3);
        //
        //        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition_4,(&crypt, zero_vector<double>(3), -unit_vector<double>(3,1)));
        //        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_4);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer_1,(&crypt, 12*crypt_radius*unit_vector<double>(3,0), unit_vector<double>(3,0)));
        simulator.AddCellKiller(p_cell_killer_1);

        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer_2,(&crypt, 12*crypt_radius*unit_vector<double>(3,1), unit_vector<double>(3,1)));
        simulator.AddCellKiller(p_cell_killer_2);

        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer_3,(&crypt, zero_vector<double>(3), -unit_vector<double>(3,0)));
        simulator.AddCellKiller(p_cell_killer_3);

        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer_4,(&crypt, zero_vector<double>(3), -unit_vector<double>(3,1)));
        simulator.AddCellKiller(p_cell_killer_4);

        // Create an instance of a Wnt concentration
        WntConcentration<3>::Instance()->SetType(LINEAR);
        WntConcentration<3>::Instance()->SetCellPopulation(crypt);
        WntConcentration<3>::Instance()->SetCryptLength(2.0*crypt_length+2.0*crypt_radius);

        // Run simulation
        simulator.Solve();

        // Label each cell according to its current node index
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();
        simulator.SetEndTime(210.0); //200

        // Run simulation
        simulator.Solve();
    }

};

#endif /*TESTCELLBASEDMULTIPLECRYPT_HPP_*/
