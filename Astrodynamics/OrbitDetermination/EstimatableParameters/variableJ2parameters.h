/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_VARIABLEJ2PARAMETERS_H
#define TUDAT_VARIABLEJ2PARAMETERS_H

#include <map>

#include <functional>

#include "Tudat/Astrodynamics/Relativity/variableJ2Interface.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{


class VariableJ2Amplitude: public EstimatableParameter< double >
{

public:

    VariableJ2Amplitude(
            const std::string& associatedBody,
            const std::shared_ptr< relativity::VariableJ2Interface > variableJ2Interface
                = relativity::variableJ2Interface):
        EstimatableParameter< double >( variable_J2_amplitude, associatedBody ),
        variableJ2Interface_( variableJ2Interface )
    { }

    //! Destructor
    ~VariableJ2Amplitude( ) { }

    double getParameterValue( )
    { return variableJ2Interface_->getAmplitude(); }

    void setParameterValue( double parameterValue )
    { variableJ2Interface_->setAmplitude(parameterValue); }

    int getParameterSize(){ return 1; }

protected:

private:

    std::shared_ptr< relativity::VariableJ2Interface > variableJ2Interface_;

};

class VariableJ2Period: public EstimatableParameter< double >
{

public:

    VariableJ2Period(
            const std::string& associatedBody,
            const std::shared_ptr< relativity::VariableJ2Interface > variableJ2Interface
                = relativity::variableJ2Interface):
        EstimatableParameter< double >( variable_J2_period, associatedBody ),
        variableJ2Interface_( variableJ2Interface )
    { }

    //! Destructor
    ~VariableJ2Period( ) { }

    double getParameterValue( )
    { return variableJ2Interface_->getPeriod(); }

    void setParameterValue( double parameterValue )
    { variableJ2Interface_->setPeriod(parameterValue); }

    int getParameterSize(){ return 1; }

protected:

private:

    std::shared_ptr< relativity::VariableJ2Interface > variableJ2Interface_;

};

class VariableJ2Phase: public EstimatableParameter< double >
{

public:

    VariableJ2Phase(
            const std::string& associatedBody,
            const std::shared_ptr< relativity::VariableJ2Interface > variableJ2Interface
                = relativity::variableJ2Interface):
        EstimatableParameter< double >( variable_J2_phase, associatedBody ),
        variableJ2Interface_( variableJ2Interface )
    { }

    //! Destructor
    ~VariableJ2Phase( ) { }

    double getParameterValue( )
    { return variableJ2Interface_->getPhase(); }

    void setParameterValue( double parameterValue )
    { variableJ2Interface_->setPhase(parameterValue); }

    int getParameterSize(){ return 1; }

protected:

private:

    std::shared_ptr< relativity::VariableJ2Interface > variableJ2Interface_;

};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_VARIABLEJ2PARAMETERS_H
