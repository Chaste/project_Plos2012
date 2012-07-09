/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "MultipleCryptGeometryBoundaryCondition.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"

MultipleCryptGeometryBoundaryCondition::MultipleCryptGeometryBoundaryCondition(AbstractCellPopulation<3>* pCellPopulation,
                                                                      double cryptRadius,
                                                                      double cryptLength,
                                                                      double villusRadius,
                                                                      double villusLength,
                                                                      double domainWidth)
    : AbstractCellPopulationBoundaryCondition<3>(pCellPopulation),
      mRadiusOfCrypt(cryptRadius),
      mLengthOfCrypt(cryptLength),
      mRadiusOfVillus(villusRadius),
      mLengthOfVillus(villusLength),
      mDomainWidth(domainWidth)
{
    assert(mRadiusOfCrypt > 0.0);

    if (dynamic_cast<NodeBasedCellPopulation<3>*>(this->mpCellPopulation) == NULL)
    {
        EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition object.");
    }
}

double MultipleCryptGeometryBoundaryCondition::GetRadiusOfCrypt() const
{
    return mRadiusOfCrypt;
}

double MultipleCryptGeometryBoundaryCondition::GetLengthOfCrypt() const
{
    return mLengthOfCrypt;
}

double MultipleCryptGeometryBoundaryCondition::GetRadiusOfVillus() const
{
    return mRadiusOfVillus;
}

double MultipleCryptGeometryBoundaryCondition::GetLengthOfVillus() const
{
    return mLengthOfVillus;
}


double MultipleCryptGeometryBoundaryCondition::GetDomainWidth() const
{
    return mDomainWidth;
}

void MultipleCryptGeometryBoundaryCondition::ImposeBoundaryCondition(const std::vector<c_vector<double, 3> >& rOldLocations)
{
    // Iterate over the cell population
    for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        const c_vector<double,3> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        //PRINT_VECTOR(cell_location);#

        const double flat_height =  2.0*mRadiusOfCrypt + mLengthOfCrypt;

        c_vector<double,3> base_centre_1;
        base_centre_1[0] = 0.5*mDomainWidth;
        base_centre_1[1] = 0.0;
        base_centre_1[2] = mRadiusOfCrypt;

        c_vector<double,3> location_from_centre_1 = cell_location - base_centre_1;
        location_from_centre_1[2]=0.0;

        c_vector<double,3> base_centre_2;
        base_centre_2[0] = 0.0;
        base_centre_2[1] = 0.5*mDomainWidth;
        base_centre_2[2] = mRadiusOfCrypt;

        c_vector<double,3> location_from_centre_2 = cell_location - base_centre_2;
        location_from_centre_2[2]=0.0;

        c_vector<double,3> base_centre_3;
        base_centre_3[0] = 0.5*mDomainWidth;
        base_centre_3[1] = mDomainWidth;
        base_centre_3[2] = mRadiusOfCrypt;

        c_vector<double,3> location_from_centre_3 = cell_location - base_centre_3;
        location_from_centre_3[2]=0.0;

        c_vector<double,3> base_centre_4;
        base_centre_4[0] = mDomainWidth;
        base_centre_4[1] = 0.5*mDomainWidth;
        base_centre_4[2] = mRadiusOfCrypt;

        c_vector<double,3> location_from_centre_4 = cell_location - base_centre_4;
        location_from_centre_4[2]=0.0;

        c_vector<double,3> base_centre_villus;
        base_centre_villus[0] = 0.5*mDomainWidth;
        base_centre_villus[1] = 0.5*mDomainWidth;
        base_centre_villus[2] = flat_height + mLengthOfVillus + mRadiusOfVillus;

        c_vector<double,3> location_from_centre_vilus = cell_location - base_centre_villus;
        location_from_centre_vilus[2]=0.0;

        c_vector<double, 3> location_on_surface = cell_location;


        // Now Crypts
        if (norm_2(location_from_centre_1) < 2.0*mRadiusOfCrypt)
        {   // 1st Crypt
            location_on_surface = MoveToCrypt(cell_location,base_centre_1);
        }
        else if (norm_2(location_from_centre_2) < 2.0*mRadiusOfCrypt)
        {   // 2nd Crypt
            location_on_surface = MoveToCrypt(cell_location,base_centre_2);
        }
        else if (norm_2(location_from_centre_3) < 2.0*mRadiusOfCrypt)
        {   // 3rd Crypt
            location_on_surface = MoveToCrypt(cell_location,base_centre_3);
        }
        else if (norm_2(location_from_centre_4) < 2.0*mRadiusOfCrypt)
        {   // 4th Crypt
            location_on_surface = MoveToCrypt(cell_location,base_centre_4);
        }
        else if (norm_2(location_from_centre_vilus) < 2.0*mRadiusOfVillus)
        {   // On villus
            location_on_surface = cell_location;
            if (cell_location[2] <= flat_height + mRadiusOfVillus) // Base of Villus
            {
                c_vector<double,3> centre_on_rim = (location_from_centre_vilus)*(2.0*mRadiusOfVillus)/norm_2(location_from_centre_vilus) + base_centre_villus;
                centre_on_rim[2] = flat_height + mRadiusOfVillus;
                double radius = norm_2(cell_location - centre_on_rim);
                location_on_surface = mRadiusOfVillus*(cell_location - centre_on_rim)/radius + centre_on_rim;
            }
            else if (   cell_location[2] > flat_height + mRadiusOfVillus
                     && cell_location[2] < flat_height + mRadiusOfVillus + mLengthOfVillus) // cylinder of villus
            {
                location_on_surface[2] = base_centre_villus[2];
                double radius = norm_2(location_on_surface - base_centre_villus);
                location_on_surface = mRadiusOfVillus*(location_on_surface - base_centre_villus)/radius + base_centre_villus;
                location_on_surface[2] = cell_location[2];
            }
            else  // Top of villus
            {
                double radius = norm_2(cell_location - base_centre_villus);
                location_on_surface = mRadiusOfVillus*(cell_location - base_centre_villus)/radius + base_centre_villus;
            }
        }
        else // Flat part
        {
            //PRINT_VECTOR(location_on_surface);
            location_on_surface[2] = flat_height;
        }

        // Edges, treat these as reflective boundaries
        if (cell_location[0] < 0.0)
        {
            location_on_surface[0] = -cell_location[0];
        }
        if (cell_location[0] > mDomainWidth)
        {
            location_on_surface[0] = mDomainWidth - (cell_location[0] - mDomainWidth);
        }
        if (cell_location[1] < 0.0)
        {
            location_on_surface[1] = -cell_location[1];
        }
        if (cell_location[1] > mDomainWidth)
        {
            location_on_surface[1] = mDomainWidth - (cell_location[1] - mDomainWidth);
        }

        // Move node on to surface
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);
        p_node->rGetModifiableLocation() = location_on_surface;
        //PRINT_VECTOR(location_on_surface);
        //PRINT_VECTOR(p_node->rGetLocation());
    }
}

c_vector<double, 3> MultipleCryptGeometryBoundaryCondition::MoveToCrypt(const c_vector<double,3>& rCellLocation, const c_vector<double,3>& rBaseCentre)
{
    c_vector<double, 3> location_on_surface = rCellLocation;
    // End of crypt
    if (rCellLocation[2] <= mRadiusOfCrypt)
    {
        double radius = norm_2(rCellLocation - rBaseCentre);
        //PRINT_VARIABLE(radius);
        location_on_surface = mRadiusOfCrypt*(rCellLocation - rBaseCentre)/radius + rBaseCentre;
    }
    else if (rCellLocation[2]  > mRadiusOfCrypt && rCellLocation[2] < mLengthOfCrypt + mRadiusOfCrypt) // cylinder of crypt
    {   // cylinder
        location_on_surface[2] = mRadiusOfCrypt;
        double radius = norm_2(location_on_surface - rBaseCentre);
        //PRINT_VARIABLE(radius);
        location_on_surface = mRadiusOfCrypt*(location_on_surface - rBaseCentre)/radius + rBaseCentre;
        location_on_surface[2] = rCellLocation[2];
    }
    else  // rim
    {
        c_vector<double, 3> location_from_centre = rCellLocation - rBaseCentre;
        location_from_centre[2]=0.0;
        c_vector<double,3> centre_on_rim = (location_from_centre)*(2.0*mRadiusOfCrypt)/norm_2(location_from_centre) + rBaseCentre;
        centre_on_rim[2] = mRadiusOfCrypt + mLengthOfCrypt;
        double radius = norm_2(rCellLocation - centre_on_rim);
        location_on_surface = mRadiusOfCrypt*(rCellLocation - centre_on_rim)/radius + centre_on_rim;
    }
    return location_on_surface;
}

bool MultipleCryptGeometryBoundaryCondition::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

//    // Iterate over the cell population
//    for (typename AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
//         cell_iter != this->mpCellPopulation->End();
//         ++cell_iter)
//    {
//        c_vector<double,3> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
//
//        c_vector<double,3> sphere_centre = original_heightero_vector<double>(3);
//        sphere_centre[2] = mRadiusOfCrypt;
//        double target_radius = mRadiusOfCrypt;
//        double original_height = cell_location[2];
//
//        if (cell_location[2] > mRadiusOfCrypt)
//        {
//            //double original_height = cell_location[2] - mRadiusOfCrypt;
//            //target_radius *= (1.0525 - 0.05*(tanh(0.5*(original_height-5.0)) - 2*tanh(1*(original_height-9.0))));
//            cell_location[2]=mRadiusOfCrypt;
//        }
//
//        // Find the radial distance between this cell and the surface
//        double radius = normnorm_2(cell_location - sphere_centre);
//
//        // If the cell is too far from the surface of the sphere...
//        if (fabs(radius - mRadiusOfCrypt) > 1e-12)
//        {
//            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
//            condition_satisfied = false;
//        }
//    }


    return condition_satisfied;
}


void MultipleCryptGeometryBoundaryCondition::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<RadiusOfCrypt>" << mRadiusOfCrypt << "</RadiusOfCrypt>\n";
    *rParamsFile << "\t\t\t<LengthOfCrypt>" << mLengthOfCrypt << "</LengthOfCrypt>\n";
    *rParamsFile << "\t\t\t<RadiusOfVillus>" << mRadiusOfVillus << "</RadiusOfVillus>\n";
    *rParamsFile << "\t\t\t<LengthOfVillus>" << mLengthOfVillus << "</LengthOfVillus>\n";
    *rParamsFile << "\t\t\t<DomainWidth>" << mDomainWidth << "</DomainWidth>\n";
    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<3>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Serialioriginal_heightation for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MultipleCryptGeometryBoundaryCondition)
