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
#ifndef TESTELECTROMECHANICSTWISTINGCUBE_HPP_
#define TESTELECTROMECHANICSTWISTINGCUBE_HPP_


#include <cxxtest/TestSuite.h>
#include "PlaneStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "CardiacElectroMechanicsProblem.hpp"

// Notes: this test requires comments. Takes 3.5 mins to run with b=GccOpt_ndebug on
// PP laptop.

class TestElectroMechanicsTwisingCube : public CxxTest::TestSuite
{
public:
    void TestTwistingCube() throw(Exception)
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-1000*1000);

        // set up two meshes of 1mm by 1mm by 1mm
        TetrahedralMesh<3,3> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01, 0.1, 0.1, 0.1);

        QuadraticMesh<3> mechanics_mesh(0.02, 0.1, 0.1, 0.1);

        // fix the nodes on Z=0
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh, 2, 0.0);

        HeartConfig::Instance()->SetSimulationDuration(50.0);

        ElectroMechanicsProblemDefinition<3> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);

//// This is how to generate a fibre file for this mesh
//        OutputFileHandler handler("FibreFile");
//        out_stream p_file = handler.OpenOutputFile("5by5by5_fibres.ortho");
//        *p_file << mechanics_mesh.GetNumElements() << "\n"; // first line is number of entries
//        for(unsigned i=0; i<mechanics_mesh.GetNumElements(); i++)
//        {
//            double X = mechanics_mesh.GetElement(i)->CalculateCentroid()(0);
//            double theta = M_PI/3 - 10*X*2*M_PI/3; // 60 degrees when X=0, -60 when X=0.1;
//            *p_file <<  "0 " << cos(theta)  << " " << sin(theta)  // first three entries are fibre direction
//                    << " 0 " << -sin(theta) << " " << cos(theta)  // next three are sheet direction
//                    << " 1 0 0\n";                                // then normal to sheet direction
//        }
//        p_file->close();

        // NOTE: the following assumes the local project folder is named 'Plot2012'
        problem_defn.SetVariableFibreSheetDirectionsFile("projects/Plos2012/test/data/5by5by5_fibres.ortho", false);

        CardiacElectroMechanicsProblem<3> problem(COMPRESSIBLE,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  &cell_factory,
                                                  &problem_defn,
                                                  "TestElectroMechanicsTwistingCube");

        problem.Solve();

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }
};

#endif /* TESTELECTROMECHANICSTWISTINGCUBE_HPP_ */
