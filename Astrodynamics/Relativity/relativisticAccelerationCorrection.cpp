/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Eigen/Geometry>
#include <iostream>

#include "Tudat/Astrodynamics/Relativity/relativisticAccelerationCorrection.h"

namespace tudat
{

namespace relativity
{

//! Function to compute a term common to several relativistic acceleration terms
double calculateRelativisticAccelerationCorrectionsCommonterm(
        const double centralBodyGravitationalParameter,
        const double relativeDistance )
{
    return centralBodyGravitationalParameter /
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT *
              relativeDistance * relativeDistance * relativeDistance );
}

//! Function to compute the Schwarzschild term of the relativistic acceleration correction.
Eigen::Vector3d calculateScharzschildGravitationalAccelerationCorrection(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& relativePosition,
        const Eigen::Vector3d& relativeVelocity,
        const double relativeDistance,
        const double commonCorrectionTerm,
        const double ppnParameterGamma,
        const double ppnParameterBeta )
{
    Eigen::Vector3d acceleration = ( 2.0 * ( ppnParameterGamma + ppnParameterBeta ) *
                      centralBodyGravitationalParameter / relativeDistance - ppnParameterGamma *
                      relativeVelocity.dot( relativeVelocity ) ) * relativePosition +
            2.0 * ( 1.0 + ppnParameterGamma ) * ( relativePosition.dot( relativeVelocity ) ) * relativeVelocity;
    return commonCorrectionTerm * acceleration;
}

//! Function to compute the Schwarzschild term of the relativistic acceleration correction.
Eigen::Vector3d calculateScharzschildGravitationalAccelerationCorrection(
        double centralBodyGravitationalParameter,
        const Eigen::Vector6d& relativeState,
        double ppnParameterGamma,
        double ppnParameterBeta )
{
    return calculateScharzschildGravitationalAccelerationCorrection(
                centralBodyGravitationalParameter, relativeState.segment( 0, 3 ),
                relativeState.segment( 3, 3 ), relativeState.segment( 0, 3 ).norm( ),
                calculateRelativisticAccelerationCorrectionsCommonterm(
                    centralBodyGravitationalParameter, relativeState.segment( 0, 3 ).norm( ) ),
                ppnParameterGamma, ppnParameterBeta );
}


Eigen::Vector3d calculateScharzschildAlphaTermsAccelerationCorrection(
        const double centralBodyGravitationalParameter,
        const double acceleratedBodyGravitationalParameter,
        const Eigen::Vector3d& relativePosition,
        const Eigen::Vector3d& relativeVelocity,
        const double relativeDistance,
        const double commonCorrectionTerm,
        const double ppnParameterAlpha1,
        const double ppnParameterAlpha2)
{
    double gravitationalParameterRatio = acceleratedBodyGravitationalParameter/centralBodyGravitationalParameter;

    Eigen::Vector3d acceleration =
            (
                (2.0 + ppnParameterAlpha1) * acceleratedBodyGravitationalParameter / relativeDistance
                - 0.5 * (6.0 + ppnParameterAlpha1 + ppnParameterAlpha2)
                    * gravitationalParameterRatio
                    * (relativeVelocity.dot(relativeVelocity))
                + 1.5 * (1.0 + ppnParameterAlpha2)
                    * gravitationalParameterRatio
                    * (relativeVelocity.dot(relativePosition/relativeDistance))
                    * (relativeVelocity.dot(relativePosition/relativeDistance))
            ) * relativePosition
            - gravitationalParameterRatio
            * (2.0 - ppnParameterAlpha1 + ppnParameterAlpha2)
            * (relativePosition.dot(relativeVelocity))
            * relativeVelocity;
    return commonCorrectionTerm * acceleration;
}



//! Function to compute the Lense-Thirring term of the relativistic acceleration correction.
Eigen::Vector3d calculateLenseThirringCorrectionAcceleration(
        const Eigen::Vector3d& relativePosition,
        const Eigen::Vector3d& relativeVelocity,
        const double relativeDistance,
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        const double ppnParameterGamma )
{
    Eigen::Vector3d termsInBrackets = 3.0
            * relativePosition.cross( relativeVelocity )
            * relativePosition.dot( centralBodyAngularMomentum )
            / (relativeDistance * relativeDistance);

    termsInBrackets += relativeVelocity.cross( centralBodyAngularMomentum );

    double commonCorrectionTerm = centralBodyGravitationalParameter
            / (physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT
               * relativeDistance * relativeDistance * relativeDistance);

    Eigen::Vector3d acceleration = ( 1.0 + ppnParameterGamma ) * commonCorrectionTerm * termsInBrackets;

//    std::cout<<"LT acc: "<<acceleration.transpose()<<std::endl;

    return acceleration;
//    return Eigen::Vector3d::Zero();
}

////! Function to compute the Lense-Thirring term of the relativistic acceleration correction.
//Eigen::Vector3d calculateLenseThirringCorrectionAcceleration(
//        const double centralBodyGravitationalParameter,
//        const Eigen::Vector6d& relativeState,
//        const Eigen::Vector3d& centralBodyAngularMomentum,
//        const double ppnParameterGamma )
//{
//    return calculateLenseThirringCorrectionAcceleration(
//                relativeState.segment( 0, 3 ),
//                relativeState.segment( 3, 3 ),
//                relativeState.segment( 0, 3 ).norm( ),
//                calculateRelativisticAccelerationCorrectionsCommonterm(
//                    centralBodyGravitationalParameter, relativeState.segment( 0, 3 ).norm( ) ),
//                centralBodyAngularMomentum, ppnParameterGamma );

//}

//! Function to compute the de Sitter term of the relativistic acceleration correction.
Eigen::Vector3d calculateDeSitterCorrectionAcceleration(
        const Eigen::Vector3d& orbiterRelativeVelocity,
        const Eigen::Vector3d& orbitedBodyPositionWrtLargerBody,
        const Eigen::Vector3d& orbitedBodyVelocityWrtLargerBody,
        const double commonCorrectionTermOfLargerBody,
        const double ppnParameterGamma )
{
    return - commonCorrectionTermOfLargerBody * ( 1.0 + 2.0 * ppnParameterGamma ) *
            ( orbitedBodyVelocityWrtLargerBody.cross( orbitedBodyPositionWrtLargerBody ) ).cross( orbiterRelativeVelocity );
}

//! Function to compute the de Sitter term of the relativistic acceleration correction.
Eigen::Vector3d calculateDeSitterCorrectionAcceleration(
        const double largerBodyGravitationalParameter,
        const Eigen::Vector6d& orbiterRelativeState,
        const Eigen::Vector6d&orbitedBodyStateWrtLargerBody,
        const double ppnParameterGamma )
{
    return calculateDeSitterCorrectionAcceleration(
                orbiterRelativeState.segment( 3, 3 ),
                orbitedBodyStateWrtLargerBody.segment( 0, 3 ),
                orbitedBodyStateWrtLargerBody.segment( 3, 3 ),
                calculateRelativisticAccelerationCorrectionsCommonterm(
                    largerBodyGravitationalParameter,
                    orbitedBodyStateWrtLargerBody.segment( 0, 3 ).norm( ) ),
                ppnParameterGamma );
}

}

}
