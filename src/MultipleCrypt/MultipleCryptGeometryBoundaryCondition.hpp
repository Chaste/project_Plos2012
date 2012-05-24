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
class MultipleCryptGeometryBoundaryCondition : public AbstractCellPopulationBoundaryCondition<3>
{
private:

    /** The radius of the crypt bases. */
    double mRadiusOfCrypt;

    /** The length of the crypts. */
    double mLengthOfCrypt;

    /** The radius of the villus base. */
    double mRadiusOfVillus;

    /** The length of the villus. */
    double mLengthOfVillus;

    /** The domain width in x and y directions */
    double mDomainWidth;

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
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<3> >(*this);
        archive & mRadiusOfCrypt;
        archive & mLengthOfCrypt;
        archive & mRadiusOfVillus;
        archive & mLengthOfVillus;
        archive & mDomainWidth;
    }

    /**
     * Move a cell onto a crypt
     *
     * @param rCellLocation  the old cell location after a movement
     * @param rBaseCentre  the location of the crypt base
     * @return the new cell location after BC enforced
     */
    c_vector<double, 3> MoveToCrypt(const c_vector<double,3>& rCellLocation, const c_vector<double,3>& rBaseCentre);


public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param cryptRadius the radius of the crypt bases (and top rim)
     * @param cryptLength the length of the crypts
     * @param villusRadius the radius of the villus top (and base rim)
     * @param villusLength the length of the villus
     * @param domainWidth the width (in x and y directions) of the domain.
     */
    MultipleCryptGeometryBoundaryCondition(AbstractCellPopulation<3>* pCellPopulation,
                                    double cryptRadius,
                                    double cryptLength,
                                    double villusRadius,
                                    double villusLength,
                                    double domainWidth);

    /**
     * @return The radii of the crypts.
     */
    double GetRadiusOfCrypt() const;

    /**
     * @return The length of the crypts.
     */
    double GetLengthOfCrypt() const;

    /**
     * @return The radius of the villus.
     */
    double GetRadiusOfVillus() const;

    /**
     * @return The length of the villus.
     */
    double GetLengthOfVillus() const;

    /**
     * @return The width of the domain (in x and y directions).
     */
    double GetDomainWidth() const;

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::vector< c_vector<double, 3> >& rOldLocations);

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
CHASTE_CLASS_EXPORT(MultipleCryptGeometryBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a MultipleCryptGeometryBoundaryCondition.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const MultipleCryptGeometryBoundaryCondition* t, const BOOST_PFTO unsigned int file_version)
{
    // Get data required to construct instance
    const AbstractCellPopulation<3>* const p_cell_population = t->GetCellPopulation();
    double crypt_radius = t->GetRadiusOfCrypt();
    double crypt_length = t->GetLengthOfCrypt();
    double villus_radius = t->GetRadiusOfVillus();
    double villus_length = t->GetLengthOfVillus();
    double domain_width = t->GetDomainWidth();

    // Archive
    ar << p_cell_population;
    ar << crypt_radius;
    ar << crypt_length;
    ar << villus_radius;
    ar << villus_length;
    ar << domain_width;
}

/**
 * De-serialize constructor parameters and initialize a MultipleCryptGeometryBoundaryCondition.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, MultipleCryptGeometryBoundaryCondition* t, const unsigned int file_version)
{
    // Make variables required to construct new instance with constructor
    AbstractCellPopulation<3>* p_cell_population;
    double crypt_radius;
    double crypt_length;
    double villus_radius;
    double villus_length;
    double domain_width;

    // Populate these variables
    ar >> p_cell_population;
    ar >> crypt_radius;
    ar >> crypt_length;
    ar >> villus_radius;
    ar >> villus_length;
    ar >> domain_width;

    // Invoke constructor to initialise instance
    ::new(t)MultipleCryptGeometryBoundaryCondition(p_cell_population, crypt_radius, crypt_length,
                                                   villus_radius, villus_length, domain_width);
}
}
} // namespace ...

#endif /*MULTIPLECRYPTGEOMETRYBOUNDARYCONDITION_HPP_*/
