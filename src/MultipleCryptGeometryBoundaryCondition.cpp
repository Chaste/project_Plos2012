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

template<unsigned DIM>
MultipleCryptGeometryBoundaryCondition<DIM>::MultipleCryptGeometryBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                                                      double radius,
                                                                      double length)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
      mRadiusOfBase(radius),
      mLengthOfCrypt(length)
{
    assert(mRadiusOfBase > 0.0);

    if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation) == NULL)
    {
        EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition object.");
    }
    if (DIM < 3)
    {
        EXCEPTION("This boundary condition is only implemented in 3D.");
    }
}

template<unsigned DIM>
double MultipleCryptGeometryBoundaryCondition<DIM>::GetRadiusOfBase() const
{
    return mRadiusOfBase;
}

template<unsigned DIM>
double MultipleCryptGeometryBoundaryCondition<DIM>::GetLengthOfCrypt() const
{
    return mLengthOfCrypt;
}

template<unsigned DIM>
void MultipleCryptGeometryBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::vector<c_vector<double, DIM> >& rOldLocations)
{
    double mLengthOfCrypt = 2.0;

    // Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        //PRINT_VECTOR(cell_location);

        double original_height = cell_location[DIM-1];

        c_vector<double,3> base_centre_1;
        base_centre_1[0] = 2.0*mRadiusOfBase;
        base_centre_1[1] = 2.0*mRadiusOfBase;
        base_centre_1[2] = mRadiusOfBase;

        c_vector<double,3> location_from_centre_1 = cell_location - base_centre_1;
        location_from_centre_1[2]=0.0;

        c_vector<double,3> base_centre_2;
        base_centre_2[0] = 10.0*mRadiusOfBase;
        base_centre_2[1] = 2.0*mRadiusOfBase;
        base_centre_2[2] = mRadiusOfBase;

        c_vector<double,3> location_from_centre_2 = cell_location - base_centre_2;
        location_from_centre_2[2]=0.0;

        c_vector<double,3> base_centre_3;
        base_centre_3[0] = 2.0*mRadiusOfBase;
        base_centre_3[1] = 10.0*mRadiusOfBase;
        base_centre_3[2] = mRadiusOfBase;

        c_vector<double,3> location_from_centre_3 = cell_location - base_centre_3;
        location_from_centre_3[2]=0.0;

        c_vector<double,3> base_centre_4;
        base_centre_4[0] = 10.0*mRadiusOfBase;
        base_centre_4[1] = 10.0*mRadiusOfBase;
        base_centre_4[2] = mRadiusOfBase;

        c_vector<double,3> location_from_centre_4 = cell_location - base_centre_4;
        location_from_centre_4[2]=0.0;

        c_vector<double,3> base_centre_vilus;
        base_centre_vilus[0] = 6.0*mRadiusOfBase;
        base_centre_vilus[1] = 6.0*mRadiusOfBase;
        base_centre_vilus[2] = 2.0*mLengthOfCrypt + 3.0* mRadiusOfBase;

        c_vector<double,3> location_from_centre_vilus = cell_location - base_centre_vilus;
        location_from_centre_vilus[2]=0.0;


        c_vector<double, DIM> location_on_surface = cell_location;

        if (norm_2(location_from_centre_1) < 2.0*mRadiusOfBase) //1st Crypt
        {
            // Base of crypt
            if (cell_location[DIM-1] <= mRadiusOfBase)
            {
                double radius = norm_2(cell_location - base_centre_1);
                //PRINT_VARIABLE(radius);
                location_on_surface = mRadiusOfBase*(cell_location - base_centre_1)/radius + base_centre_1;
            }
            else if (original_height > mRadiusOfBase && original_height < mLengthOfCrypt + mRadiusOfBase) // cylinder of crypt
            {
                cell_location[DIM-1] = mRadiusOfBase;
                double radius = norm_2(cell_location - base_centre_1);
                //PRINT_VARIABLE(radius);
                location_on_surface = mRadiusOfBase*(cell_location - base_centre_1)/radius + base_centre_1;
                location_on_surface[DIM-1] = original_height;
            }
            else  // Top rim
            {
                c_vector<double,DIM> centre_on_rim = (location_from_centre_1)*(2.0*mRadiusOfBase)/norm_2(location_from_centre_1) + base_centre_1;
                centre_on_rim[DIM-1] = mRadiusOfBase + mLengthOfCrypt;
                double radius = norm_2(cell_location - centre_on_rim);
                location_on_surface = mRadiusOfBase*(cell_location - centre_on_rim)/radius + centre_on_rim;
            }
        }

        else if (norm_2(location_from_centre_2) < 2.0*mRadiusOfBase) //2nd Crypt
        {
            // Base of crypt
            if (cell_location[DIM-1] <= mRadiusOfBase)
            {
                double radius = norm_2(cell_location - base_centre_2);
                location_on_surface = mRadiusOfBase*(cell_location - base_centre_2)/radius + base_centre_2;
            }
            else if (original_height > mRadiusOfBase && original_height < mLengthOfCrypt + mRadiusOfBase) // cylinder of crypt
            {
                cell_location[DIM-1] = mRadiusOfBase;
                double radius = norm_2(cell_location - base_centre_2);
                location_on_surface = mRadiusOfBase*(cell_location - base_centre_2)/radius + base_centre_2;
                location_on_surface[DIM-1] = original_height;
            }
            else  // Top rim
            {
                c_vector<double,DIM> centre_on_rim = (location_from_centre_2)*(2.0*mRadiusOfBase)/norm_2(location_from_centre_2) + base_centre_2;
                centre_on_rim[DIM-1] = mRadiusOfBase + mLengthOfCrypt;
                double radius = norm_2(cell_location - centre_on_rim);
                location_on_surface = mRadiusOfBase*(cell_location - centre_on_rim)/radius + centre_on_rim;
            }
        }

        else if (norm_2(location_from_centre_3) < 2.0*mRadiusOfBase) //3rd Crypt
        {
            // Base of crypt
            if (cell_location[DIM-1] <= mRadiusOfBase)
            {
                double radius = norm_2(cell_location - base_centre_3);
                location_on_surface = mRadiusOfBase*(cell_location - base_centre_3)/radius + base_centre_3;
            }
            else if (original_height > mRadiusOfBase && original_height < mLengthOfCrypt + mRadiusOfBase) // cylinder of crypt
            {
                cell_location[DIM-1] = mRadiusOfBase;
                double radius = norm_2(cell_location - base_centre_3);
                location_on_surface = mRadiusOfBase*(cell_location - base_centre_3)/radius + base_centre_3;
                location_on_surface[DIM-1] = original_height;
            }
            else  // Top rim
            {
                c_vector<double,DIM> centre_on_rim = (location_from_centre_3)*(2.0*mRadiusOfBase)/norm_2(location_from_centre_3) + base_centre_3;
                centre_on_rim[DIM-1] = mRadiusOfBase + mLengthOfCrypt;
                double radius = norm_2(cell_location - centre_on_rim);
                location_on_surface = mRadiusOfBase*(cell_location - centre_on_rim)/radius + centre_on_rim;
            }
        }

        else if (norm_2(location_from_centre_4) < 2.0*mRadiusOfBase) //4th Crypt
        {
            // Base of crypt
            if (cell_location[DIM-1] <= mRadiusOfBase)
            {
                double radius = norm_2(cell_location - base_centre_4);
                location_on_surface = mRadiusOfBase*(cell_location - base_centre_4)/radius + base_centre_4;
            }
            else if (original_height > mRadiusOfBase && original_height < mLengthOfCrypt + mRadiusOfBase) // cylinder of crypt
            {
                cell_location[DIM-1] = mRadiusOfBase;
                double radius = norm_2(cell_location - base_centre_4);
                location_on_surface = mRadiusOfBase*(cell_location - base_centre_4)/radius + base_centre_4;
                location_on_surface[DIM-1] = original_height;
            }
            else  // Top rim
            {
                c_vector<double,DIM> centre_on_rim = (location_from_centre_4)*(2.0*mRadiusOfBase)/norm_2(location_from_centre_4) + base_centre_4;
                centre_on_rim[DIM-1] = mRadiusOfBase + mLengthOfCrypt;
                double radius = norm_2(cell_location - centre_on_rim);
                location_on_surface = mRadiusOfBase*(cell_location - centre_on_rim)/radius + centre_on_rim;
            }
        }

        else if (norm_2(location_from_centre_vilus) < 2.0*mRadiusOfBase) //Vilus
        {
            // Top of Vilus
            if (cell_location[DIM-1] >= base_centre_vilus[DIM-1])
            {
                double radius = norm_2(cell_location - base_centre_vilus);
                location_on_surface = mRadiusOfBase*(cell_location - base_centre_vilus)/radius + base_centre_vilus;
            }
            else if (original_height < base_centre_vilus[DIM-1] && original_height > base_centre_vilus[DIM-1] - mLengthOfCrypt) // cylinder of vilus                        {
            {
                cell_location[DIM-1] = base_centre_vilus[DIM-1];
                double radius = norm_2(cell_location - base_centre_vilus);
                location_on_surface = mRadiusOfBase*(cell_location - base_centre_vilus)/radius + base_centre_vilus;
                location_on_surface[DIM-1] = original_height;
            }
            else  // Top rim
            {
                c_vector<double,DIM> centre_on_rim = (location_from_centre_vilus)*(2.0*mRadiusOfBase)/norm_2(location_from_centre_vilus) + base_centre_vilus;
                centre_on_rim[DIM-1] = 3.0*mRadiusOfBase + mLengthOfCrypt;
                double radius = norm_2(cell_location - centre_on_rim);
                location_on_surface = mRadiusOfBase*(cell_location - centre_on_rim)/radius + centre_on_rim;
            }
        }


        else // Flat part
        {
            //PRINT_VECTOR(location_on_surface);
            location_on_surface[DIM-1] =  2.0*mRadiusOfBase + mLengthOfCrypt;
        }

        // Move node on to surface
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
        p_node->rGetModifiableLocation() = location_on_surface;
        //PRINT_VECTOR(location_on_surface);
        //PRINT_VECTOR(p_node->rGetLocation());

    }


}

template<unsigned DIM>
bool MultipleCryptGeometryBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

//    // Iterate over the cell population
//    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
//         cell_iter != this->mpCellPopulation->End();
//         ++cell_iter)
//    {
//        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
//
//        c_vector<double,DIM> sphere_centre = original_heightero_vector<double>(DIM);
//        sphere_centre[DIM-1] = mRadiusOfBase;
//        double target_radius = mRadiusOfBase;
//        double original_height = cell_location[DIM-1];
//
//        if (cell_location[DIM-1] > mRadiusOfBase)
//        {
//            //double original_height = cell_location[DIM-1] - mRadiusOfBase;
//            //target_radius *= (1.0525 - 0.05*(tanh(0.5*(original_height-5.0)) - 2*tanh(1*(original_height-9.0))));
//            cell_location[DIM-1]=mRadiusOfBase;
//        }
//
//        // Find the radial distance between this cell and the surface
//        double radius = normnorm_2(cell_location - sphere_centre);
//
//        // If the cell is too far from the surface of the sphere...
//        if (fabs(radius - mRadiusOfBase) > 1e-12)
//        {
//            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
//            condition_satisfied = false;
//        }
//    }


    return condition_satisfied;
}

template<unsigned DIM>
void MultipleCryptGeometryBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<RadiusOfBase>" << mRadiusOfBase << "</RadiusOfBase>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MultipleCryptGeometryBoundaryCondition<1>;
template class MultipleCryptGeometryBoundaryCondition<2>;
template class MultipleCryptGeometryBoundaryCondition<3>;

// Serialioriginal_heightation for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MultipleCryptGeometryBoundaryCondition)
