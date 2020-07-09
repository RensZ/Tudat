/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_VARIABLEJ4PARAMETERS_H
#define TUDAT_VARIABLEJ4PARAMETERS_H

#include <map>

#include <functional>

#include "Tudat/Astrodynamics/Relativity/variableJ4Interface.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{


class VariableJ4Amplitude: public EstimatableParameter< double >
{

public:

    VariableJ4Amplitude(
            const std::string& associatedBody,
            const std::shared_ptr< relativity::VariableJ4Interface > variableJ4Interface
                = relativity::variableJ4Interface):
        EstimatableParameter< double >( variable_J4_amplitude, associatedBody ),
        variableJ4Interface_( variableJ4Interface )
    { }

    //! Destructor
    ~VariableJ4Amplitude( ) { }

    double getParameterValue( )
    { return variableJ4Interface_->getAmplitude(); }

    void setParameterValue( double parameterValue )
    { variableJ4Interface_->setAmplitude(parameterValue); }

    int getParameterSize(){ return 1; }

protected:

private:

    std::shared_ptr< relativity::VariableJ4Interface > variableJ4Interface_;

};

class VariableJ4Period: public EstimatableParameter< double >
{

public:

    VariableJ4Period(
            const std::string& associatedBody,
            const std::shared_ptr< relativity::VariableJ4Interface > variableJ4Interface
                = relativity::variableJ4Interface):
        EstimatableParameter< double >( variable_J4_period, associatedBody ),
        variableJ4Interface_( variableJ4Interface )
    { }

    //! Destructor
    ~VariableJ4Period( ) { }

    double getParameterValue( )
    { return variableJ4Interface_->getPeriod(); }

    void setParameterValue( double parameterValue )
    { variableJ4Interface_->setPeriod(parameterValue); }

    int getParameterSize(){ return 1; }

protected:

private:

    std::shared_ptr< relativity::VariableJ4Interface > variableJ4Interface_;

};

class VariableJ4Phase: public EstimatableParameter< double >
{

public:

    VariableJ4Phase(
            const std::string& associatedBody,
            const std::shared_ptr< relativity::VariableJ4Interface > variableJ4Interface
                = relativity::variableJ4Interface):
        EstimatableParameter< double >( variable_J4_phase, associatedBody ),
        variableJ4Interface_( variableJ4Interface )
    { }

    //! Destructor
    ~VariableJ4Phase( ) { }

    double getParameterValue( )
    { return variableJ4Interface_->getPhase(); }

    void setParameterValue( double parameterValue )
    { variableJ4Interface_->setPhase(parameterValue); }

    int getParameterSize(){ return 1; }

protected:

private:

    std::shared_ptr< relativity::VariableJ4Interface > variableJ4Interface_;

};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_VARIABLEJ4PARAMETERS_H
