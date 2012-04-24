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

#ifndef MULTIPLECRYPTGEOMETRYBOUNDARYCONDITION_HPP_
#define MULTIPLECRYPTGEOMETRYBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A cell population boundary condition class, which restricts nodes to lie
 * on the surface of a set of crypts and a villus.
 */
template<unsigned DIM>
class MultipleCryptGeometryBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /** The radius of the base. */
    double mRadiusOfBase;

    /** The length of the crypt. */
    double mLengthOfCrypt;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param radius the radius of the base
     * @param length the length of the crypt
     */
    MultipleCryptGeometryBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                    double radius,
                                    double length);

    /**
     * @return mRadiusOfBase.
     */
    double GetRadiusOfBase() const;

    /**
     * @return mLengthOfCrypt.
     */
    double GetLengthOfCrypt() const;

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::vector< c_vector<double, DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MultipleCryptGeometryBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a MultipleCryptGeometryBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MultipleCryptGeometryBoundaryCondition<DIM>* t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    // Archive other member variables
    double radius = t->GetRadiusOfBase();
    ar << radius;

    double length = t->GetLengthOfCrypt();
    ar << length;
}

/**
 * De-serialize constructor parameters and initialize a MultipleCryptGeometryBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MultipleCryptGeometryBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Retrieve other member variables
    double radius;
    ar >> radius;
    double length;
    ar >> length;

    // Invoke inplace constructor to initialise instance
    ::new(t)MultipleCryptGeometryBoundaryCondition<DIM>(p_cell_population, radius, length);
}
}
} // namespace ...

#endif /*MULTIPLECRYPTGEOMETRYBOUNDARYCONDITION_HPP_*/
