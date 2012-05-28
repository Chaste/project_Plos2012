#ifndef TESTCELLBASEDMULTIPLECRYPT_HPP_
#define TESTCELLBASEDMULTIPLECRYPT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "CellsGenerator.hpp"
#include "RandomCellKiller.hpp"
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
        double crypt_length = 4.0;
        double crypt_radius = 1.0;
        double villus_length = 8.0;
        double villus_radius = 2.0;
        double domain_width = 12.0;
        double domain_height = 2.0*crypt_radius+crypt_length+villus_length+2.0*villus_radius;

        // Put a single cell at the base of each crypt.
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u,  false,  0.25*domain_width, 0.25*domain_width, 0.0));
        nodes.push_back(new Node<3>(2u,  false,  0.25*domain_width, 0.75*domain_width, 0.0));
        nodes.push_back(new Node<3>(4u,  false,  0.75*domain_width, 0.25*domain_width, 0.0));
        nodes.push_back(new Node<3>(6u,  false,  0.75*domain_width, 0.75*domain_width, 0.0));

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
        simulator.SetOutputDirectory("Plos2012_MultipleCrypt");
        simulator.SetDt(1.0/180.0);
        simulator.SetSamplingTimestepMultiple(60);
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

        MAKE_PTR_ARGS(RandomCellKiller<3>, p_cell_killer_2,(&crypt, 0.008)); // prob of death in an hour
        simulator.AddCellKiller(p_cell_killer_2);

        // Create an instance of a Wnt concentration
        WntConcentration<3>::Instance()->SetType(LINEAR);
        WntConcentration<3>::Instance()->SetCellPopulation(crypt);
        WntConcentration<3>::Instance()->SetCryptLength(crypt_length+2*crypt_radius);

        // Run simulation
        simulator.SetEndTime(200);
        simulator.Solve(); // to 200

        // Label each cell according to its current node index
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();
        simulator.SetEndTime(1000.0);

        // Run simulation
        simulator.Solve();
    }

};

#endif /*TESTCELLBASEDMULTIPLECRYPT_HPP_*/
