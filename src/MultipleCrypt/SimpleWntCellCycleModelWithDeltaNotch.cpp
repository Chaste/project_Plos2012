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

#include "UblasIncludes.hpp"
#include "SimpleWntCellCycleModelWithDeltaNotch.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "Exception.hpp"


SimpleWntCellCycleModelWithDeltaNotch::SimpleWntCellCycleModelWithDeltaNotch(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : CellCycleModelOdeHandler(DOUBLE_UNSET, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<SimpleWntCellCycleModelWithDeltaNotch, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<SimpleWntCellCycleModelWithDeltaNotch, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

AbstractCellCycleModel* SimpleWntCellCycleModelWithDeltaNotch::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    SimpleWntCellCycleModelWithDeltaNotch* p_model = new SimpleWntCellCycleModelWithDeltaNotch(this->mpOdeSolver);

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide) will already have been
     * correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetDimension(mDimension);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetUseCellProliferativeTypeDependentG1Duration(mUseCellProliferativeTypeDependentG1Duration);
    p_model->SetWntStemThreshold(mWntStemThreshold);
    p_model->SetWntTransitThreshold(mWntTransitThreshold);
    p_model->SetWntLabelledThreshold(mWntLabelledThreshold);
    p_model->SetWntLabelledThreshold(mWntLabelledThreshold);
    p_model->SetLastTime(mLastTime);

    // Create the new cell-cycle model's ODE system
    p_model->SetOdeSystem(new DeltaNotchOdeSystem());

    // Use the current values of the state variables in mpOdeSystem as an initial condition for the new cell-cycle model's ODE system
    assert(mpOdeSystem);
    p_model->SetStateVariables(mpOdeSystem->rGetStateVariables());


    return p_model;
}

void SimpleWntCellCycleModelWithDeltaNotch::UpdateCellCyclePhase()
{
    assert(SimulationTime::Instance()->IsStartTimeSetUp());
    UpdateDeltaNotch();
    SolveOdeToTime(SimulationTime::Instance()->GetTime());
    SimpleWntCellCycleModel::UpdateCellCyclePhase();

    this->mpCell->GetCellData()->SetItem("notch", GetNotch());
    this->mpCell->GetCellData()->SetItem("delta", GetDelta());
}

void SimpleWntCellCycleModelWithDeltaNotch::Initialise()
{
    assert(mpOdeSystem == NULL);
    assert(mpCell != NULL);

    mpOdeSystem = new DeltaNotchOdeSystem;
    if(mInitialConditions == std::vector<double>())
    {
        mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
    }
    else
    {
        mpOdeSystem->SetStateVariables(mInitialConditions);
    }

    SimpleWntCellCycleModel::Initialise();

    SetLastTime(mBirthTime);
}


void SimpleWntCellCycleModelWithDeltaNotch::SetInitialConditions(std::vector<double> initialConditions)
{
    assert(initialConditions.size() == 2);
    mInitialConditions = initialConditions;
}

void SimpleWntCellCycleModelWithDeltaNotch::UpdateDeltaNotch()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);

    double mean_delta = mpCell->GetCellData()->GetItem("mean delta");
    mpOdeSystem->SetParameter("Mean Delta", mean_delta);
}

double SimpleWntCellCycleModelWithDeltaNotch::GetNotch()
{
    assert(mpOdeSystem != NULL);
    double notch = mpOdeSystem->rGetStateVariables()[0];
    return notch;
}

double SimpleWntCellCycleModelWithDeltaNotch::GetDelta()
{
    assert(mpOdeSystem != NULL);
    double delta = mpOdeSystem->rGetStateVariables()[1];
    return delta;
}

double SimpleWntCellCycleModelWithDeltaNotch::GetMeanNeighbouringDelta()
{
    assert(mpOdeSystem != NULL);
    double mean_neighbouring_delta = mpOdeSystem->GetParameter("Mean Delta");
    return mean_neighbouring_delta;
}

void SimpleWntCellCycleModelWithDeltaNotch::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output.

    // Call direct parent class
    SimpleWntCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleWntCellCycleModelWithDeltaNotch)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(SimpleWntCellCycleModelWithDeltaNotch)
