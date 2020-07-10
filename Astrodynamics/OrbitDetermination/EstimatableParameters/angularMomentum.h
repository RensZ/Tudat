/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ANGULARMOMENTUM_H
#define TUDAT_ANGULARMOMENTUM_H

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a angular momentum
class AngularMomentum: public EstimatableParameter< double >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param rotationModel rotation model object containing the angular momentum to be estimated.
     * \param associatedBody Name of body containing the gravityFieldModel object
     */
    AngularMomentum(
            const std::shared_ptr< ephemerides::RotationalEphemeris > rotationModel, const std::string& associatedBody ):
        EstimatableParameter< double >( angular_momentum, associatedBody ),
        rotationModel_( rotationModel ){ }

    //! Destructor
    ~AngularMomentum( ) { }

    //! Function to get the current value of the angular momentum that is to be estimated.
    /*!
     * Function to get the current value of the angular momentum that is to be estimated.
     * \return Current value of the angular momentum that is to be estimated.
     */
    double getParameterValue( )
    {
        return rotationModel_->getAngularMomentum( );
    }

    //! Function to reset the value of the angular momentum that is to be estimated.
    /*!
     * Function to reset the value of the angular momentum that is to be estimated.
     * \param parameterValue New value of the angular momentum that is to be estimated.
     */
    void setParameterValue( double parameterValue )
    {
        rotationModel_->resetAngularMomentum( parameterValue );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    //! rotation model object containing the angular momentum to be estimated.
    std::shared_ptr< ephemerides::RotationalEphemeris > rotationModel_;

};

} // namespace estimatable_parameters

} // namespace tudat


#endif // TUDAT_GRAVITATIONALPARAMETER_H
