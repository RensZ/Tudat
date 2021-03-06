/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>
#include <functional>
#include <memory>

#include <boost/make_shared.hpp>
#include <boost/bind.hpp>
#include "Tudat/Astrodynamics/Aerodynamics/flightConditions.h"
#include "Tudat/Astrodynamics/Ephemerides/frameManager.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/Propulsion/thrustMagnitudeWrapper.h"
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Astrodynamics/Relativity/relativisticAccelerationCorrection.h"
#include "Tudat/Astrodynamics/Relativity/metric.h"
#include "Tudat/Basics/utilities.h"
#include "Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createFlightConditions.h"

#include "tudatApplications/thesis/MyApplications/timeVaryingGravitationalParameterAcceleration.h"
#include "tudatApplications/thesis/MyApplications/sepViolationAcceleration.h"

#include "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/customFunctions.h"

//#include "tudatApplications/thesis/MyApplications/timeVaryingGravitationalParameter.h"
//#include "tudatApplications/thesis/MyApplications/TVGPInterface.h"
//#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"

namespace tudat
{

namespace simulation_setup
{

using namespace aerodynamics;
using namespace gravitation;
using namespace basic_astrodynamics;
using namespace electro_magnetism;
using namespace ephemerides;


//! Function to create a direct (i.e. not third-body) gravitational acceleration (of any type)
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > createDirectGravitationalAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfCentralBody,
        const bool isCentralBody )
{
    // Check if sum of gravitational parameters (i.e. inertial force w.r.t. central body) should be used.
    bool sumGravitationalParameters = 0;
    if( ( nameOfCentralBody == nameOfBodyExertingAcceleration ) && bodyUndergoingAcceleration != nullptr )
    {
        sumGravitationalParameters = 1;
    }


    // Check type of acceleration model and create.
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel;
    switch( accelerationSettings->accelerationType_ )
    {
    case central_gravity:
        accelerationModel = createCentralGravityAcceleratioModel(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    sumGravitationalParameters );
        break;
    case spherical_harmonic_gravity:
        accelerationModel = createSphericalHarmonicsGravityAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings,
                    sumGravitationalParameters );
        break;
    case mutual_spherical_harmonic_gravity:
        accelerationModel = createMutualSphericalHarmonicsGravityAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings,
                    sumGravitationalParameters,
                    isCentralBody );
        break;
    default:

        std::string errorMessage = "Error when making gravitional acceleration model, cannot parse type " +
                std::to_string( accelerationSettings->accelerationType_ );
        throw std::runtime_error( errorMessage );
    }
    return accelerationModel;
}

//! Function to create a third-body gravitational acceleration (of any type)
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > createThirdBodyGravitationalAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Check type of acceleration model and create.
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel;
    switch( accelerationSettings->accelerationType_ )
    {
    case central_gravity:
        accelerationModel = std::make_shared< ThirdBodyCentralGravityAcceleration >(
                    std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                        createDirectGravitationalAcceleration(
                            bodyUndergoingAcceleration, bodyExertingAcceleration,
                            nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 0 ) ),
                    std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                        createDirectGravitationalAcceleration(
                            centralBody, bodyExertingAcceleration,
                            nameOfCentralBody, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 1 ) ), nameOfCentralBody );
        break;
    case spherical_harmonic_gravity:
        accelerationModel = std::make_shared< ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                    std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            bodyUndergoingAcceleration, bodyExertingAcceleration,
                            nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 0 ) ),
                    std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            centralBody, bodyExertingAcceleration, nameOfCentralBody, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 1 ) ), nameOfCentralBody );
        break;
    case mutual_spherical_harmonic_gravity:
        accelerationModel = std::make_shared< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >(
                    std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            bodyUndergoingAcceleration, bodyExertingAcceleration,
                            nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 0 ) ),
                    std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                        createDirectGravitationalAcceleration(
                            centralBody, bodyExertingAcceleration, nameOfCentralBody, nameOfBodyExertingAcceleration,
                            accelerationSettings, "", 1 ) ), nameOfCentralBody );
        break;
    default:

        std::string errorMessage = "Error when making third-body gravitional acceleration model, cannot parse type " +
                std::to_string( accelerationSettings->accelerationType_ );
        throw std::runtime_error( errorMessage );
    }
    return accelerationModel;
}

//! Function to create gravitational acceleration (of any type)
std::shared_ptr< AccelerationModel< Eigen::Vector3d > > createGravitationalAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfCentralBody )
{

    std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelPointer;
    if( accelerationSettings->accelerationType_ != central_gravity &&
            accelerationSettings->accelerationType_ != spherical_harmonic_gravity &&
            accelerationSettings->accelerationType_ != mutual_spherical_harmonic_gravity )
    {
        throw std::runtime_error( "Error when making gravitational acceleration, type is inconsistent" );
    }

    if( nameOfCentralBody == nameOfBodyExertingAcceleration || ephemerides::isFrameInertial( nameOfCentralBody ) )
    {
        accelerationModelPointer = createDirectGravitationalAcceleration( bodyUndergoingAcceleration,
                                                                          bodyExertingAcceleration,
                                                                          nameOfBodyUndergoingAcceleration,
                                                                          nameOfBodyExertingAcceleration,
                                                                          accelerationSettings,
                                                                          nameOfCentralBody, false );
    }
    else
    {
        accelerationModelPointer = createThirdBodyGravitationalAcceleration( bodyUndergoingAcceleration,
                                                                             bodyExertingAcceleration,
                                                                             centralBody,
                                                                             nameOfBodyUndergoingAcceleration,
                                                                             nameOfBodyExertingAcceleration,
                                                                             nameOfCentralBody, accelerationSettings );
    }

    return accelerationModelPointer;
}


//! Function to create central gravity acceleration model.
std::shared_ptr< CentralGravitationalAccelerationModel3d > createCentralGravityAcceleratioModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const bool useCentralBodyFixedFrame )
{
    // Declare pointer to return object.
    std::shared_ptr< CentralGravitationalAccelerationModel3d > accelerationModelPointer;

    // Check if body is endowed with a gravity field model (i.e. is capable of exerting
    // gravitation acceleration).
    if( bodyExertingAcceleration->getGravityFieldModel( ) == nullptr )
    {
        throw std::runtime_error(
                    std::string( "Error, gravity field model not set when making central ") +
                    " gravitational acceleration of " + nameOfBodyExertingAcceleration + " on " +
                    nameOfBodyUndergoingAcceleration );
    }
    else
    {
        std::function< double( ) > gravitationalParameterFunction;

        // Set correct value for gravitational parameter.
        if( useCentralBodyFixedFrame == 0  ||
                bodyUndergoingAcceleration->getGravityFieldModel( ) == nullptr )
        {
            gravitationalParameterFunction =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodyExertingAcceleration->getGravityFieldModel( ) );
        }
        else
        {
            std::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodyExertingAcceleration->getGravityFieldModel( ) );
            std::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodyUndergoingAcceleration->getGravityFieldModel( ) );
            gravitationalParameterFunction =
                    std::bind( &utilities::sumFunctionReturn< double >,
                               gravitationalParameterOfBodyExertingAcceleration,
                               gravitationalParameterOfBodyUndergoingAcceleration );
        }

        // Create acceleration object.
        accelerationModelPointer =
                std::make_shared< CentralGravitationalAccelerationModel3d >(
                    std::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                    gravitationalParameterFunction,
                    std::bind( &Body::getPosition, bodyExertingAcceleration ),
                    useCentralBodyFixedFrame );
    }


    return accelerationModelPointer;
}

//! Function to create spherical harmonic gravity acceleration model.
std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel >
createSphericalHarmonicsGravityAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame,
        const bool useDegreeZeroTerm )
{
    // Declare pointer to return object
    std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicsSettings =
            std::dynamic_pointer_cast< SphericalHarmonicAccelerationSettings >(
                accelerationSettings );
    if( sphericalHarmonicsSettings == nullptr )
    {
        throw std::runtime_error(
                    std::string( "Error, acceleration settings inconsistent ") +
                    " making sh gravitational acceleration of " + nameOfBodyExertingAcceleration +
                    " on " + nameOfBodyUndergoingAcceleration );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );

        std::shared_ptr< RotationalEphemeris> rotationalEphemeris =
                bodyExertingAcceleration->getRotationalEphemeris( );
        if( sphericalHarmonicsGravityField == nullptr )
        {
            throw std::runtime_error(
                        std::string( "Error, spherical harmonic gravity field model not set when ")
                        + " making sh gravitational acceleration of " +
                        nameOfBodyExertingAcceleration +
                        " on " + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            if( rotationalEphemeris == nullptr )
            {
                throw std::runtime_error( "Warning when making spherical harmonic acceleration on body " +
                                          nameOfBodyUndergoingAcceleration + ", no rotation model found for " +
                                          nameOfBodyExertingAcceleration );
            }

            if( rotationalEphemeris->getTargetFrameOrientation( ) !=
                    sphericalHarmonicsGravityField->getFixedReferenceFrame( ) )
            {
                throw std::runtime_error( "Warning when making spherical harmonic acceleration on body " +
                                          nameOfBodyUndergoingAcceleration + ", rotation model found for " +
                                          nameOfBodyExertingAcceleration + " is incompatible, frames are: " +
                                          rotationalEphemeris->getTargetFrameOrientation( ) + " and " +
                                          sphericalHarmonicsGravityField->getFixedReferenceFrame( ) );
            }

            std::function< double( ) > gravitationalParameterFunction;

            // Check if mutual acceleration is to be used.
            if( useCentralBodyFixedFrame == false ||
                    bodyUndergoingAcceleration->getGravityFieldModel( ) == nullptr )
            {
                gravitationalParameterFunction =
                        std::bind( &SphericalHarmonicsGravityField::getGravitationalParameter,
                                   sphericalHarmonicsGravityField );
            }
            else
            {
                // Create function returning summed gravitational parameter of the two bodies.
                std::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                   sphericalHarmonicsGravityField );
                std::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                   bodyUndergoingAcceleration->getGravityFieldModel( ) );
                gravitationalParameterFunction =
                        std::bind( &utilities::sumFunctionReturn< double >,
                                   gravitationalParameterOfBodyExertingAcceleration,
                                   gravitationalParameterOfBodyUndergoingAcceleration );
            }

            std::function< Eigen::MatrixXd( ) > originalCosineCoefficientFunction =
                    std::bind( &SphericalHarmonicsGravityField::getCosineCoefficientsBlock,
                               sphericalHarmonicsGravityField,
                               sphericalHarmonicsSettings->maximumDegree_,
                               sphericalHarmonicsSettings->maximumOrder_ );

            std::function< Eigen::MatrixXd( ) > cosineCoefficientFunction;
            if( !useDegreeZeroTerm )
            {
                cosineCoefficientFunction =
                        std::bind( &setDegreeAndOrderCoefficientToZero, originalCosineCoefficientFunction );
            }
            else
            {
                cosineCoefficientFunction = originalCosineCoefficientFunction;
            }

            // Create acceleration object.
            accelerationModel =
                    std::make_shared< SphericalHarmonicsGravitationalAccelerationModel >
                    ( std::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                      gravitationalParameterFunction,
                      sphericalHarmonicsGravityField->getReferenceRadius( ),
                      cosineCoefficientFunction,
                      std::bind( &SphericalHarmonicsGravityField::getSineCoefficientsBlock,
                                 sphericalHarmonicsGravityField,
                                 sphericalHarmonicsSettings->maximumDegree_,
                                 sphericalHarmonicsSettings->maximumOrder_ ),
                      std::bind( &Body::getPosition, bodyExertingAcceleration ),
                      std::bind( &Body::getCurrentRotationToGlobalFrame,
                                 bodyExertingAcceleration ), useCentralBodyFixedFrame );
        }
    }
    return accelerationModel;
}

//! Function to create mutual spherical harmonic gravity acceleration model.
std::shared_ptr< gravitation::MutualSphericalHarmonicsGravitationalAccelerationModel >
createMutualSphericalHarmonicsGravityAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame,
        const bool acceleratedBodyIsCentralBody )
{
    using namespace basic_astrodynamics;

    // Declare pointer to return object
    std::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< MutualSphericalHarmonicAccelerationSettings > mutualSphericalHarmonicsSettings =
            std::dynamic_pointer_cast< MutualSphericalHarmonicAccelerationSettings >( accelerationSettings );
    if( mutualSphericalHarmonicsSettings == nullptr )
    {
        std::string errorMessage = "Error, expected mutual spherical harmonics acceleration settings when making acceleration model on " +
                nameOfBodyUndergoingAcceleration + "due to " + nameOfBodyExertingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyExertingAcceleration =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyUndergoingAcceleration->getGravityFieldModel( ) );

        if( sphericalHarmonicsGravityFieldOfBodyExertingAcceleration == nullptr )
        {

            std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );

        }
        else if( sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration == nullptr )
        {

            std::string errorMessage = "Error " + nameOfBodyUndergoingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::function< double( ) > gravitationalParameterFunction;

            // Create function returning summed gravitational parameter of the two bodies.
            if( useCentralBodyFixedFrame == false )
            {
                gravitationalParameterFunction =
                        std::bind( &SphericalHarmonicsGravityField::getGravitationalParameter,
                                   sphericalHarmonicsGravityFieldOfBodyExertingAcceleration );
            }
            else
            {
                // Create function returning summed gravitational parameter of the two bodies.
                std::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                   sphericalHarmonicsGravityFieldOfBodyExertingAcceleration );
                std::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                   sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration );
                gravitationalParameterFunction =
                        std::bind( &utilities::sumFunctionReturn< double >,
                                   gravitationalParameterOfBodyExertingAcceleration,
                                   gravitationalParameterOfBodyUndergoingAcceleration );
            }

            // Create acceleration object.

            int maximumDegreeOfUndergoingBody, maximumOrderOfUndergoingBody;
            if( !acceleratedBodyIsCentralBody )
            {
                maximumDegreeOfUndergoingBody = mutualSphericalHarmonicsSettings->maximumDegreeOfBodyUndergoingAcceleration_;
                maximumOrderOfUndergoingBody = mutualSphericalHarmonicsSettings->maximumOrderOfBodyUndergoingAcceleration_;
            }
            else
            {
                maximumDegreeOfUndergoingBody = mutualSphericalHarmonicsSettings->maximumDegreeOfCentralBody_;
                maximumOrderOfUndergoingBody = mutualSphericalHarmonicsSettings->maximumOrderOfCentralBody_;
            }

            accelerationModel = std::make_shared< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                        std::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                        std::bind( &Body::getPosition, bodyExertingAcceleration ),
                        gravitationalParameterFunction,
                        sphericalHarmonicsGravityFieldOfBodyExertingAcceleration->getReferenceRadius( ),
                        sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration->getReferenceRadius( ),
                        std::bind( &SphericalHarmonicsGravityField::getCosineCoefficientsBlock,
                                   sphericalHarmonicsGravityFieldOfBodyExertingAcceleration,
                                   mutualSphericalHarmonicsSettings->maximumDegreeOfBodyExertingAcceleration_,
                                   mutualSphericalHarmonicsSettings->maximumOrderOfBodyExertingAcceleration_ ),
                        std::bind( &SphericalHarmonicsGravityField::getSineCoefficientsBlock,
                                   sphericalHarmonicsGravityFieldOfBodyExertingAcceleration,
                                   mutualSphericalHarmonicsSettings->maximumDegreeOfBodyExertingAcceleration_,
                                   mutualSphericalHarmonicsSettings->maximumOrderOfBodyExertingAcceleration_ ),
                        std::bind( &SphericalHarmonicsGravityField::getCosineCoefficientsBlock,
                                   sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration,
                                   maximumDegreeOfUndergoingBody,
                                   maximumOrderOfUndergoingBody ),
                        std::bind( &SphericalHarmonicsGravityField::getSineCoefficientsBlock,
                                   sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration,
                                   maximumDegreeOfUndergoingBody,
                                   maximumOrderOfUndergoingBody ),
                        std::bind( &Body::getCurrentRotationToGlobalFrame,
                                   bodyExertingAcceleration ),
                        std::bind( &Body::getCurrentRotationToGlobalFrame,
                                   bodyUndergoingAcceleration ),
                        useCentralBodyFixedFrame );
        }
    }
    return accelerationModel;
}


//! Function to create a third body central gravity acceleration model.
std::shared_ptr< gravitation::ThirdBodyCentralGravityAcceleration >
createThirdBodyCentralGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody )
{
    // Declare pointer to return object.
    std::shared_ptr< ThirdBodyCentralGravityAcceleration > accelerationModelPointer;

    // Create acceleration object.
    accelerationModelPointer =  std::make_shared< ThirdBodyCentralGravityAcceleration >(
                std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                    createCentralGravityAcceleratioModel( bodyUndergoingAcceleration,
                                                          bodyExertingAcceleration,
                                                          nameOfBodyUndergoingAcceleration,
                                                          nameOfBodyExertingAcceleration, 0 ) ),
                std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                    createCentralGravityAcceleratioModel( centralBody, bodyExertingAcceleration,
                                                          nameOfCentralBody,
                                                          nameOfBodyExertingAcceleration, 0 ) ), nameOfCentralBody );

    return accelerationModelPointer;
}

//! Function to create a third body spheric harmonic gravity acceleration model.
std::shared_ptr< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >
createThirdBodySphericalHarmonicGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    using namespace basic_astrodynamics;

    // Declare pointer to return object
    std::shared_ptr< ThirdBodySphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicsSettings =
            std::dynamic_pointer_cast< SphericalHarmonicAccelerationSettings >( accelerationSettings );
    if( sphericalHarmonicsSettings == nullptr )
    {
        std::string errorMessage = "Error, expected spherical harmonics acceleration settings when making acceleration model on " +
                nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );
        if( sphericalHarmonicsGravityField == nullptr )
        {
            std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making third body spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else
        {

            accelerationModel =  std::make_shared< ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                        std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                            createSphericalHarmonicsGravityAcceleration(
                                bodyUndergoingAcceleration, bodyExertingAcceleration, nameOfBodyUndergoingAcceleration,
                                nameOfBodyExertingAcceleration, sphericalHarmonicsSettings, 0 ) ),
                        std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                            createSphericalHarmonicsGravityAcceleration(
                                centralBody, bodyExertingAcceleration, nameOfCentralBody,
                                nameOfBodyExertingAcceleration, sphericalHarmonicsSettings, 0 ) ), nameOfCentralBody );
        }
    }
    return accelerationModel;
}

//! Function to create a third body mutual spheric harmonic gravity acceleration model.
std::shared_ptr< gravitation::ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >
createThirdBodyMutualSphericalHarmonicGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Declare pointer to return object
    std::shared_ptr< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< MutualSphericalHarmonicAccelerationSettings > mutualSphericalHarmonicsSettings =
            std::dynamic_pointer_cast< MutualSphericalHarmonicAccelerationSettings >( accelerationSettings );
    if( mutualSphericalHarmonicsSettings == nullptr )
    {

        std::string errorMessage = "Error, expected mutual spherical harmonics acceleration settings when making acceleration model on " +
                nameOfBodyUndergoingAcceleration +
                " due to " + nameOfBodyExertingAcceleration;
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyExertingAcceleration =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyUndergoingAcceleration->getGravityFieldModel( ) );
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityFieldOfCentralBody =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    centralBody->getGravityFieldModel( ) );

        if( sphericalHarmonicsGravityFieldOfBodyExertingAcceleration == nullptr )
        {
            std::string errorMessage = "Error " + nameOfBodyExertingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else if( sphericalHarmonicsGravityFieldOfBodyUndergoingAcceleration == nullptr )
        {
            std::string errorMessage = "Error " + nameOfBodyUndergoingAcceleration + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else if( sphericalHarmonicsGravityFieldOfCentralBody == nullptr )
        {
            std::string errorMessage = "Error " + nameOfCentralBody + " does not have a spherical harmonics gravity field " +
                    "when making mutual spherical harmonics gravity acceleration on " +
                    nameOfBodyUndergoingAcceleration;
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::shared_ptr< MutualSphericalHarmonicAccelerationSettings > accelerationSettingsForCentralBodyAcceleration =
                    std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                        mutualSphericalHarmonicsSettings->maximumDegreeOfBodyExertingAcceleration_,
                        mutualSphericalHarmonicsSettings->maximumOrderOfBodyExertingAcceleration_,
                        mutualSphericalHarmonicsSettings->maximumDegreeOfCentralBody_,
                        mutualSphericalHarmonicsSettings->maximumOrderOfCentralBody_ );
            accelerationModel =  std::make_shared< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >(
                        std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                            createMutualSphericalHarmonicsGravityAcceleration(
                                bodyUndergoingAcceleration, bodyExertingAcceleration, nameOfBodyUndergoingAcceleration,
                                nameOfBodyExertingAcceleration, mutualSphericalHarmonicsSettings, 0, 0 ) ),
                        std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                            createMutualSphericalHarmonicsGravityAcceleration(
                                centralBody, bodyExertingAcceleration, nameOfCentralBody,
                                nameOfBodyExertingAcceleration, accelerationSettingsForCentralBodyAcceleration, 0, 1 ) ),
                        nameOfCentralBody );
        }
    }
    return accelerationModel;
}

//! Function to create an aerodynamic acceleration model.
std::shared_ptr< aerodynamics::AerodynamicAcceleration > createAerodynamicAcceleratioModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    // Check existence of required environment models
    if( bodyUndergoingAcceleration->getAerodynamicCoefficientInterface( ) == nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic acceleration, body " +
                                  nameOfBodyUndergoingAcceleration +
                                  "has no aerodynamic coefficients." );
    }

    if( bodyExertingAcceleration->getAtmosphereModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic acceleration, central body " +
                                  nameOfBodyExertingAcceleration + " has no atmosphere model.");
    }

    if( bodyExertingAcceleration->getShapeModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic acceleration, central body " +
                                  nameOfBodyExertingAcceleration + " has no shape model." );
    }

    // Retrieve flight conditions; create object if not yet extant.
    std::shared_ptr< AtmosphericFlightConditions > bodyFlightConditions =
            std::dynamic_pointer_cast< AtmosphericFlightConditions >( bodyUndergoingAcceleration->getFlightConditions( ) );

    if( bodyFlightConditions == nullptr && bodyUndergoingAcceleration->getFlightConditions( ) == nullptr )
    {
        bodyFlightConditions = createAtmosphericFlightConditions( bodyUndergoingAcceleration,
                                                                  bodyExertingAcceleration,
                                                                  nameOfBodyUndergoingAcceleration,
                                                                  nameOfBodyExertingAcceleration );
        bodyUndergoingAcceleration->setFlightConditions( bodyFlightConditions );
    }
    else if( bodyFlightConditions == nullptr && bodyUndergoingAcceleration->getFlightConditions( ) != nullptr )
    {
        throw std::runtime_error( "Error when making aerodynamic acceleration, found flight conditions that are not atmospheric." );
    }

    // Retrieve frame in which aerodynamic coefficients are defined.
    std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > aerodynamicCoefficients =
            bodyUndergoingAcceleration->getAerodynamicCoefficientInterface( );
    reference_frames::AerodynamicsReferenceFrames accelerationFrame;
    if( aerodynamicCoefficients->getAreCoefficientsInAerodynamicFrame( ) )
    {
        accelerationFrame = reference_frames::aerodynamic_frame;
    }
    else
    {
        accelerationFrame = reference_frames::body_frame;
    }

    // Create function to transform from frame of aerodynamic coefficienrs to that of propagation.
    std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > toPropagationFrameTransformation;
    toPropagationFrameTransformation =
            reference_frames::getAerodynamicForceTransformationFunction(
                bodyFlightConditions->getAerodynamicAngleCalculator( ),
                accelerationFrame,
                std::bind( &Body::getCurrentRotationToGlobalFrame, bodyExertingAcceleration ),
                reference_frames::inertial_frame );

    std::function< Eigen::Vector3d( ) > coefficientFunction =
            std::bind( &AerodynamicCoefficientInterface::getCurrentForceCoefficients,
                       aerodynamicCoefficients );
    std::function< Eigen::Vector3d( ) > coefficientInPropagationFrameFunction =
            std::bind( &reference_frames::transformVectorFunctionFromVectorFunctions,
                       coefficientFunction, toPropagationFrameTransformation );

    // Create acceleration model.
    return std::make_shared< AerodynamicAcceleration >(
                coefficientInPropagationFrameFunction,
                std::bind( &AtmosphericFlightConditions::getCurrentDensity, bodyFlightConditions ),
                std::bind( &AtmosphericFlightConditions::getCurrentAirspeed, bodyFlightConditions ),
                std::bind( &Body::getBodyMass, bodyUndergoingAcceleration ),
                std::bind( &AerodynamicCoefficientInterface::getReferenceArea,
                           aerodynamicCoefficients ),
                aerodynamicCoefficients->getAreCoefficientsInNegativeAxisDirection( ) );
}

//! Function to create a cannonball radiation pressure acceleration model.
std::shared_ptr< CannonBallRadiationPressureAcceleration >
createCannonballRadiationPressureAcceleratioModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    // Retrieve radiation pressure interface
    if( bodyUndergoingAcceleration->getRadiationPressureInterfaces( ).count(
                nameOfBodyExertingAcceleration ) == 0 )
    {
        throw std::runtime_error(
                    "Error when making radiation pressure, no radiation pressure interface found  in " +
                    nameOfBodyUndergoingAcceleration +
                    " for body " + nameOfBodyExertingAcceleration );
    }
    std::shared_ptr< RadiationPressureInterface > radiationPressureInterface =
            bodyUndergoingAcceleration->getRadiationPressureInterfaces( ).at(
                nameOfBodyExertingAcceleration );

    // Create acceleration model.
    return std::make_shared< CannonBallRadiationPressureAcceleration >(
                std::bind( &Body::getPosition, bodyExertingAcceleration ),
                std::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                std::bind( &RadiationPressureInterface::getCurrentRadiationPressure, radiationPressureInterface ),
                std::bind( &RadiationPressureInterface::getRadiationPressureCoefficient, radiationPressureInterface ),
                std::bind( &RadiationPressureInterface::getArea, radiationPressureInterface ),
                std::bind( &Body::getBodyMass, bodyUndergoingAcceleration ) );

}

//! Function to create a panelled radiation pressure acceleration model.
std::shared_ptr< electro_magnetism::PanelledRadiationPressureAcceleration > createPanelledRadiationPressureAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration )
{
    using namespace tudat::electro_magnetism;

    // Declare pointer to return object.
    std::shared_ptr< PanelledRadiationPressureAcceleration > accelerationModel;

    // Get radiation pressure interface from body undergoing acceleration, containing data on how body responds to radiation pressure.
    std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface =
            std::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                bodyUndergoingAcceleration->getRadiationPressureInterfaces( ).at( nameOfBodyExertingAcceleration ) );

    if( radiationPressureInterface == NULL )
    {
        throw std::runtime_error(
                    "Error, body undergoing acceleration, " + nameOfBodyUndergoingAcceleration +
                    " possesses no radiation pressure coefficient interface when making panelled radiation pressure acceleration due to " +
                    nameOfBodyExertingAcceleration );
    }
    else
    {
        // Create acceleration model.
        accelerationModel = std::make_shared< PanelledRadiationPressureAcceleration >(
                    radiationPressureInterface, std::bind( &Body::getBodyMass, bodyUndergoingAcceleration ) );
    }
    return accelerationModel;
}

//! Function to create a solar sail radiation pressure acceleration model.
std::shared_ptr< SolarSailAcceleration > createSolarSailAccelerationModel(
    const std::shared_ptr< Body > bodyUndergoingAcceleration,
    const std::shared_ptr< Body > bodyExertingAcceleration,
    const std::shared_ptr< Body > centralBody,
    const std::string& nameOfBodyUndergoingAcceleration,
    const std::string& nameOfBodyExertingAcceleration )
{
    // Retrieve radiation pressure interface.
    if( bodyUndergoingAcceleration->getRadiationPressureInterfaces( ).count(
            nameOfBodyExertingAcceleration ) == 0 )
    {
        throw std::runtime_error(
            "Error when making radiation pressure, no radiation pressure interface found  in " +
            nameOfBodyUndergoingAcceleration +
            " for body " + nameOfBodyExertingAcceleration );
    }

    // Get radiation pressure interface from body undergoing acceleration, containing data on how body responds to radiation pressure.
    std::shared_ptr< SolarSailingRadiationPressureInterface > radiationPressureInterface =
            std::dynamic_pointer_cast< SolarSailingRadiationPressureInterface >(
                bodyUndergoingAcceleration->getRadiationPressureInterfaces( ).at( nameOfBodyExertingAcceleration ) );

    // Create and return solar sailing acceleration model.
    return std::make_shared< SolarSailAcceleration >(
                std::bind( &Body::getPosition, bodyExertingAcceleration ),
                std::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                std::bind( &Body::getVelocity, bodyUndergoingAcceleration ),
                std::bind( &Body::getVelocity, centralBody ),
                std::bind( &SolarSailingRadiationPressureInterface::getCurrentRadiationPressure, radiationPressureInterface ),
                std::bind( &SolarSailingRadiationPressureInterface::getCurrentConeAngle, radiationPressureInterface ),
                std::bind( &SolarSailingRadiationPressureInterface::getCurrentClockAngle, radiationPressureInterface ),
                std::bind( &SolarSailingRadiationPressureInterface::getFrontEmissivityCoefficient, radiationPressureInterface ),
                std::bind( &SolarSailingRadiationPressureInterface::getBackEmissivityCoefficient, radiationPressureInterface ),
                std::bind( &SolarSailingRadiationPressureInterface::getFrontLambertianCoefficient, radiationPressureInterface ),
                std::bind( &SolarSailingRadiationPressureInterface::getBackLambertianCoefficient, radiationPressureInterface ),
                std::bind( &SolarSailingRadiationPressureInterface::getReflectivityCoefficient, radiationPressureInterface ),
                std::bind( &SolarSailingRadiationPressureInterface::getSpecularReflectionCoefficient, radiationPressureInterface ),
                std::bind( &RadiationPressureInterface::getArea, radiationPressureInterface ),
                std::bind( &Body::getBodyMass, bodyUndergoingAcceleration ) );

}


//! Function to calculate acceleration due to a time-varying graviational parameter
std::shared_ptr< TimeVaryingGravitationalParameterAcceleration > createTimeVaryingGravitationalParameterAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings)
{
    using namespace relativity;

    if (nameOfBodyExertingAcceleration != "Sun"){
        throw std::runtime_error( "Warning, time varying gravitational parameter is not tested for other bodies than the Sun, current body exerting acceleration: "
                                  + nameOfBodyExertingAcceleration );
    }

    // Declare pointer to return object
    std::shared_ptr< TimeVaryingGravitationalParameterAcceleration > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< TimeVaryingGravitationalParameterAccelerationSettings > timeVaryingGravitationalParameterAccelerationSettings =
            std::dynamic_pointer_cast< TimeVaryingGravitationalParameterAccelerationSettings >(
                accelerationSettings );
    if( timeVaryingGravitationalParameterAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error, expected TVGP acceleration settings when making acceleration model on " +
                                  nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration );
    }
    else
    {


//        std::shared_ptr< TVGPInterface > timeVarying =
//                bodyUndergoingAcceleration->getAerodynamicCoefficientInterface( );

        // Retrieve function pointers for properties of bodies exerting/undergoing acceleration.
        std::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingAcceleration =
                std::bind( &Body::getState, bodyExertingAcceleration );
        std::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingAcceleration =
                std::bind( &Body::getState, bodyUndergoingAcceleration );

        std::function< double( ) > centralBodyGravitationalParameterFunction;
        centralBodyGravitationalParameterFunction =
                std::bind( &GravityFieldModel::getGravitationalParameter, bodyExertingAcceleration->getGravityFieldModel( ) );

        ppnParameterSet->setTimeVaryingGravitationalParameter(
                    timeVaryingGravitationalParameterAccelerationSettings->timeVaryingGravitationalParameter_);

        std::function< double( ) >  timeVaryingGravitationalParameterFunction;
        timeVaryingGravitationalParameterFunction =
                std::bind( &PPNParameterSet::getTimeVaryingGravitationalParameter, ppnParameterSet );

        accelerationModel = std::make_shared< TimeVaryingGravitationalParameterAcceleration >
                ( stateFunctionOfBodyUndergoingAcceleration,
                  stateFunctionOfBodyExertingAcceleration,
                  centralBodyGravitationalParameterFunction,
                  timeVaryingGravitationalParameterFunction
                  );

    }
    return accelerationModel;
}




//! get the gravitational self energy of a body, which is needed in the function below
//! values and equation come from Genova et al. 2018, nature communications
//! radii from //https://nssdc.gsfc.nasa.gov/planetary/factsheet/<INSERTPLANET>fact.html
double getGravitationalSelfEnergy(
        const std::string bodyName,
        const double gravitationalParameter)
{
    double gravitationalSelfEnergy;

    double mass = gravitationalParameter/physical_constants::GRAVITATIONAL_CONSTANT;

    // bodies for which gravitational self energy value is quite well known
    if (bodyName == "Sun"){
        gravitationalSelfEnergy = -3.52E-6
                *mass
                *physical_constants::SPEED_OF_LIGHT
                *physical_constants::SPEED_OF_LIGHT;
    }
    else if (bodyName == "Earth"){
        gravitationalSelfEnergy = -4.64E-10
                *mass
                *physical_constants::SPEED_OF_LIGHT
                *physical_constants::SPEED_OF_LIGHT;
    }
    else if (bodyName == "Moon"){
        gravitationalSelfEnergy = -1.88E-11
                *mass
                *physical_constants::SPEED_OF_LIGHT
                *physical_constants::SPEED_OF_LIGHT;
    }

    // for the other planets, approximate by assuming a sphere with uniform density
    else{

        double bodyRadius;
        if (bodyName == "Mercury"){ bodyRadius = 2439.7E3; }
        else if (bodyName == "Venus"){ bodyRadius = 6051.8E3; }
        else if (bodyName == "Mars"){ bodyRadius = 3389.5E3; }
        else if (bodyName == "Jupiter"){ bodyRadius = 69911.0E3; }
        else if (bodyName == "Saturn"){ bodyRadius = 58232.0E3; }
        else if (bodyName == "Uranus"){ bodyRadius = 25362.0E3; }
        else if (bodyName == "Neptune"){ bodyRadius = 24622.0E3; }
        else if (bodyName == "2000001"){ bodyRadius = 939.4E3/2.0; } //Ceres
        else if (bodyName == "2000002"){ bodyRadius = 512.0E3/2.0; } //Pallas
        else if (bodyName == "2000004"){ bodyRadius = 525.4E3/2.0; } //Vesta
        else if (bodyName == "2000010"){ bodyRadius = 434.0E3/2.0; } //Hygiea
        else {
            throw std::runtime_error( "Error, gravitational self energy for body "
                                      + bodyName + " can not be implemented, body radius unknown");
        }

        gravitationalSelfEnergy =
                (-3.0/5.0) *
                (gravitationalParameter*mass)
                /bodyRadius;
    }

    return gravitationalSelfEnergy;
}


//! Function to get the new position of the Sun due to SEP effects
//!     This function is called when creating the acceleration due to the SEP in the function below
//!     The equation is given in Genova et al. 2018, nature communication, equation (9)
//!     A simplified result has been derived and the original equations are below
Eigen::Vector3d getSEPCorrectedPosition(
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const NamedBodyMap& bodyMap,
        std::vector< std::string > bodyNames,
        std::function< double( ) > nordtvedtParameterFunction)
{
    using namespace relativity;

    // Simplified solution by using taylor series expansion
    Eigen::Vector3d SEPPositionCorrectionSimplified;

    double nordtvedtParameter = nordtvedtParameterFunction();

    if (nordtvedtParameter == 0.0){ // to prevent dividing by 0 in the calculations below
        SEPPositionCorrectionSimplified << 0.0, 0.0, 0.0;
    }
    else {

        // get term dependent on properties of the Sun
        const double gravitationalParameterBodyExertingAcceleration =
                bodyExertingAcceleration->getGravityFieldModel( )->getGravitationalParameter();

        const double gravitationalSelfEnergyBodyExertingAcceleration =
                getGravitationalSelfEnergy(nameOfBodyExertingAcceleration,
                                           gravitationalParameterBodyExertingAcceleration);



        Eigen::Vector3d simplifiedSummationTerm = Eigen::Vector3d::Zero();

        const int numberOfBodies = bodyNames.size();
        for( int i = 0; i < numberOfBodies; i++ ){

           std::string currentBodyName = bodyNames.at(i);

           if (currentBodyName != nameOfBodyExertingAcceleration){

               const std::shared_ptr<Body> currentBody = bodyMap.at(currentBodyName);

               Eigen::Vector3d currentBodyPosition = currentBody->getPosition();

               double currentBodyGravitationalParameter =
                       currentBody->getGravityFieldModel()->getGravitationalParameter();

               double currentBodyGravitationalSelfEnergy =
                       getGravitationalSelfEnergy(currentBodyName, currentBodyGravitationalParameter);

               Eigen::Vector3d currentSimplifiedBodyTerm =
                       currentBodyGravitationalParameter
                       * currentBodyPosition
                       * ( gravitationalSelfEnergyBodyExertingAcceleration
                         /( gravitationalParameterBodyExertingAcceleration/physical_constants::GRAVITATIONAL_CONSTANT )
                         -
                         currentBodyGravitationalSelfEnergy
                         /( currentBodyGravitationalParameter/physical_constants::GRAVITATIONAL_CONSTANT )
                       );

               simplifiedSummationTerm += currentSimplifiedBodyTerm;

           }

        }

        SEPPositionCorrectionSimplified =
                (-1.0/gravitationalParameterBodyExertingAcceleration)
                * nordtvedtParameter
                * (1.0/(physical_constants::SPEED_OF_LIGHT*physical_constants::SPEED_OF_LIGHT))
                * simplifiedSummationTerm;
//        std::cout<<"mu sun: "<<gravitationalParameterBodyExertingAcceleration;
//        std::cout<<" eta: "<<nordtvedtParameter;
//        std::cout<<std::setprecision(15)<<" dr_SEP: "<<SEPPositionCorrectionSimplified.transpose()<<std::endl;
    }

//    std::cout<<SEPPositionCorrectionSimplified.transpose()<<std::endl;

//    std::cout<<"retreived SEP Sun position correction: "<<SEPPositionCorrectionSimplified.transpose()<<std::endl;

    return SEPPositionCorrectionSimplified;

}


//! Function to get the new position of the Sun due to SEP effects
//!     This function is called when creating the acceleration due to the SEP in the function below
//!     The equation is given in Genova et al. 2018, nature communication, equation (9)
//!     This is the unsimplified version, only used for verification purposes
Eigen::Vector3d getSEPCorrectedPositionUnsimplified(
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const NamedBodyMap& bodyMap,
        std::vector< std::string > bodyNames,
        std::function< double( ) > nordtvedtParameterFunction)
{

    using namespace relativity;

    // get SEP corrected position
    Eigen::Vector3d SEPPositionCorrection;

    double nordtvedtParameter = nordtvedtParameterFunction();

    if (nordtvedtParameter == 0.0){ // to prevent dividing by 0 in the calculations below
        SEPPositionCorrection << 0.0, 0.0, 0.0;
    }
    else {

        // get term dependent on properties of the Sun
        std::function< double( ) > centralBodyGravitationalParameterFunction;
        std::shared_ptr< GravityFieldModel > gravityFieldCentralBody = bodyExertingAcceleration->getGravityFieldModel( );

        const double gravitationalParameterBodyExertingAcceleration =
                gravityFieldCentralBody->getGravitationalParameter();

        const double gravitationalSelfEnergyBodyExertingAcceleration =
                getGravitationalSelfEnergy(nameOfBodyExertingAcceleration,
                                           gravitationalParameterBodyExertingAcceleration);

        double referenceCentralBodyTerm =
                -1.0/gravitationalParameterBodyExertingAcceleration;

        double centralBodyTerm =
                -1.0/(
                    gravitationalParameterBodyExertingAcceleration *
                    (1.0 - nordtvedtParameter *
                         gravitationalSelfEnergyBodyExertingAcceleration /
                         ( ( gravitationalParameterBodyExertingAcceleration/physical_constants::GRAVITATIONAL_CONSTANT)
                           * physical_constants::SPEED_OF_LIGHT*physical_constants::SPEED_OF_LIGHT
                         )
                    )
                );

        // get term dependent on properties of the other bodies
        Eigen::Vector3d summationTerm = Eigen::Vector3d::Zero();
        Eigen::Vector3d referenceSummationTerm = Eigen::Vector3d::Zero();

        const int numberOfBodies = bodyNames.size();

        for( int i = 0; i < numberOfBodies; i++ ){

           std::string currentBodyName = bodyNames.at(i);
           Eigen::Vector3d currentBodyTerm;
           Eigen::Vector3d currentReferenceBodyTerm;

           if (currentBodyName != nameOfBodyExertingAcceleration){

               const std::shared_ptr<Body> currentBody = bodyMap.at(currentBodyName);

               Eigen::Vector3d currentBodyPosition = currentBody->getPosition();

               std::shared_ptr< gravitation::GravityFieldModel> currentBodyGravityField =
                       currentBody->getGravityFieldModel();
               double currentBodyGravitationalParameter =
                       currentBodyGravityField->getGravitationalParameter();

               double currentBodyGravitationalSelfEnergy =
                       getGravitationalSelfEnergy(currentBodyName, currentBodyGravitationalParameter);

               currentReferenceBodyTerm =
                       currentBodyGravitationalParameter
                       * currentBodyPosition;

               currentBodyTerm =
                       (1.0 - nordtvedtParameter*
                            currentBodyGravitationalSelfEnergy/
                            ( ( currentBodyGravitationalParameter/physical_constants::GRAVITATIONAL_CONSTANT)
                              * physical_constants::SPEED_OF_LIGHT*physical_constants::SPEED_OF_LIGHT
                            )
                       )
                       * currentBodyGravitationalParameter
                       * currentBodyPosition;

               referenceSummationTerm += currentReferenceBodyTerm;
               summationTerm += currentBodyTerm;

           }

        }

        // multiply first and second term to get the SEP corrected position
        SEPPositionCorrection = centralBodyTerm*summationTerm -
                referenceCentralBodyTerm*referenceSummationTerm;

    }

    return SEPPositionCorrection;
}




//! the partial equation of the SEP violation acceleration w.r.t. the Nordtvedt parameter is calculated
//!     This is passed on via sepViolationAcceleration to sepViolationAccelerationPartial
//!     And is calculated using Genova et al. 2018, nature communication, equation (10)
Eigen::Vector3d getNordtvedtPartial(
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const NamedBodyMap& bodyMap,
        std::vector< std::string > bodyNames)
{
    using namespace relativity;

    Eigen::Vector3d term1, term2, term3;
    term1 = Eigen::Vector3d::Zero();
    term2 = Eigen::Vector3d::Zero();
    term3 = Eigen::Vector3d::Zero();

    // Properties of body exerting acceleration
    Eigen::Vector3d positionBodyExertingAcceleration = bodyExertingAcceleration->getPosition( );

    const double gravitationalParameterBodyExertingAcceleration =
            bodyExertingAcceleration->getGravityFieldModel( )->getGravitationalParameter();

    const double gravitationalSelfEnergyBodyExertingAcceleration =
            getGravitationalSelfEnergy(nameOfBodyExertingAcceleration, gravitationalParameterBodyExertingAcceleration);

    // Properties of body undergoing acceleration
    Eigen::Vector3d positionBodyUndergoingAcceleration = bodyUndergoingAcceleration->getPosition( );

    Eigen::Vector3d positionBodyUndergoingAccelerationWrtBodyExertingAcceleration =
            positionBodyUndergoingAcceleration - positionBodyExertingAcceleration;
    double distanceBodyUndergoingAccelerationWrtBodyExertingAcceleration =
            positionBodyUndergoingAccelerationWrtBodyExertingAcceleration.norm();

    double gravitationalParameterBodyUndergoingAcceleration =
            bodyUndergoingAcceleration->getGravityFieldModel()->getGravitationalParameter( );

    double gravitationalSelfEnergyBodyUndergoingAcceleration =
            getGravitationalSelfEnergy(nameOfBodyUndergoingAcceleration,
                                       gravitationalParameterBodyUndergoingAcceleration);

    // Loop over bodies
    const int numberOfBodies = bodyNames.size();
    for( int i = 0; i < numberOfBodies; i++ ){

        // current body properties
        std::string currentBodyName = bodyNames.at(i);
        const std::shared_ptr<Body> currentBody = bodyMap.at(currentBodyName);

        Eigen::Vector3d positionCurrentBody = currentBody->getPosition( );

        Eigen::Vector3d positionBodyUndergoingAccelerationWrtCurrentBody =
                positionBodyUndergoingAcceleration - positionCurrentBody;
        double distanceBodyUndergoingAccelerationWrtCurrentBody =
                positionBodyUndergoingAccelerationWrtCurrentBody.norm();

        Eigen::Vector3d positionBodyExertingAccelerationWrtCurrentBody =
                positionBodyExertingAcceleration - positionCurrentBody;
        double distanceBodyExertingAccelerationWrtCurrentBody =
                positionBodyExertingAccelerationWrtCurrentBody.norm();

        double gravitationalParameterCurrentBody =
                currentBody->getGravityFieldModel()->getGravitationalParameter( );

        double gravitationalSelfEnergyCurrentBody =
                getGravitationalSelfEnergy(currentBodyName,
                                           gravitationalParameterCurrentBody);

        // first term
        if (currentBodyName != nameOfBodyUndergoingAcceleration){
            term1 -= gravitationalParameterCurrentBody
                    * gravitationalSelfEnergyBodyUndergoingAcceleration
                    * 1.0/(gravitationalParameterBodyUndergoingAcceleration/physical_constants::GRAVITATIONAL_CONSTANT)
                    * 1.0/(physical_constants::SPEED_OF_LIGHT*physical_constants::SPEED_OF_LIGHT)
                    * positionBodyUndergoingAccelerationWrtCurrentBody
                    / (distanceBodyUndergoingAccelerationWrtCurrentBody
                       * distanceBodyUndergoingAccelerationWrtCurrentBody
                       * distanceBodyUndergoingAccelerationWrtCurrentBody);
        }

        if (currentBodyName != nameOfBodyExertingAcceleration){

            // second term
            term2 -= gravitationalParameterCurrentBody
                    * gravitationalSelfEnergyBodyExertingAcceleration
                    * 1.0/(gravitationalParameterBodyExertingAcceleration/physical_constants::GRAVITATIONAL_CONSTANT)
                    * 1.0/(physical_constants::SPEED_OF_LIGHT*physical_constants::SPEED_OF_LIGHT)
                    * positionBodyExertingAccelerationWrtCurrentBody
                    / (distanceBodyExertingAccelerationWrtCurrentBody
                       * distanceBodyExertingAccelerationWrtCurrentBody
                       * distanceBodyExertingAccelerationWrtCurrentBody);

            // third term
            term3 -= gravitationalParameterCurrentBody
                    * (gravitationalSelfEnergyCurrentBody
                       * 1.0/(gravitationalParameterCurrentBody/physical_constants::GRAVITATIONAL_CONSTANT)
                       * 1.0/(physical_constants::SPEED_OF_LIGHT*physical_constants::SPEED_OF_LIGHT)
                        -
                       gravitationalSelfEnergyBodyExertingAcceleration
                       * 1.0/(gravitationalParameterBodyExertingAcceleration/physical_constants::GRAVITATIONAL_CONSTANT)
                       * 1.0/(physical_constants::SPEED_OF_LIGHT*physical_constants::SPEED_OF_LIGHT)
                       )
                    * (1.0 / (distanceBodyUndergoingAccelerationWrtBodyExertingAcceleration
                              * distanceBodyUndergoingAccelerationWrtBodyExertingAcceleration
                              * distanceBodyUndergoingAccelerationWrtBodyExertingAcceleration)
                       )
                    * (positionBodyUndergoingAccelerationWrtBodyExertingAcceleration
                       * positionBodyUndergoingAccelerationWrtBodyExertingAcceleration.transpose()
                       / (distanceBodyUndergoingAccelerationWrtBodyExertingAcceleration
                          * distanceBodyUndergoingAccelerationWrtBodyExertingAcceleration)
                       - Eigen::Matrix3d::Identity()
                       )
                    * positionCurrentBody;

        }

    }

    Eigen::Vector3d NordtvedtPartial = term1 + term2 + term3;

//    std::cout<<"retreived Nordtvedt Partial: "<<NordtvedtPartial.transpose()<<std::endl;

//    std::cout<<NordtvedtPartial.transpose()<<" = "
//             <<term1.transpose()<<" + "
//             <<term2.transpose()<<" + "
//             <<term3.transpose()<<std::endl;

    return NordtvedtPartial;

//        // Milani's approach (Milani et al 2002, Physical Review D, derivative of eq A5 w.r.t. eta)
//        if (currentBodyName != nameOfBodyExertingAcceleration && currentBodyName != nameOfBodyUndergoingAcceleration){

//            Eigen::Vector3d currentBodyPositionWrtBodyExertingAcceleration =
//                    currentBodyPosition - bodyExertingAccelerationPosition;
//            double currentBodyDistanceWrtBodyExertingAcceleration =
//                    currentBodyPositionWrtBodyExertingAcceleration.norm( );

//            NordtvedtPartial += currentBodyGravitationalParameter
//                    * currentBodyPositionWrtBodyExertingAcceleration
//                    / (currentBodyDistanceWrtBodyExertingAcceleration
//                       * currentBodyDistanceWrtBodyExertingAcceleration
//                       * currentBodyDistanceWrtBodyExertingAcceleration);
//        }

//        NordtvedtPartial *= gravitationalSelfEnergyBodyExertingAcceleration;

}



//! Function to create an acceleration correction because of a SEP violation (nonzero nordtvedt parameter)
std::shared_ptr< relativity::SEPViolationAcceleration > createSEPViolationAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const NamedBodyMap& bodyMap )
{
    using namespace relativity;

    // Declare pointer to return object
    std::shared_ptr< SEPViolationAcceleration > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< SEPViolationAccelerationSettings > sepViolationAccelerationSettings =
            std::dynamic_pointer_cast< SEPViolationAccelerationSettings >(
                accelerationSettings );
    if( sepViolationAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error, no SEP acceleration model settings given when making acceleration model on " +
                                  nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration );
    }
    else
    {


        // Retrieve function pointers for properties of bodies exerting/undergoing acceleration.
        std::function< Eigen::Vector3d( ) > positionFunctionOfBodyExertingAcceleration =
                std::bind( &Body::getPosition, bodyExertingAcceleration );
        std::function< Eigen::Vector3d( ) > positionFunctionOfBodyUndergoingAcceleration =
                std::bind( &Body::getPosition, bodyUndergoingAcceleration );


        // ppn parameters
        std::function< double( ) >  nordtvedtParameterFunction;
        if (sepViolationAccelerationSettings->useNordtvedtConstraint_
                && sepViolationAccelerationSettings->ignoreNordtvedtConstraintInEstimation_){
            nordtvedtParameterFunction =
                    std::bind( &PPNParameterSet::getNordtvedtParameterFromPpnParameters, ppnParameterSet );
        }
        else{
            nordtvedtParameterFunction =
                    std::bind( &PPNParameterSet::getNordtvedtParameter, ppnParameterSet );
        }

        std::function< bool( ) > useNordtvedtConstraintFunction;
        useNordtvedtConstraintFunction = [ = ]( ){ return
                    sepViolationAccelerationSettings->useNordtvedtConstraint_; };


        // Gravitational parameter central body
        std::function< double( ) > centralBodyGravitationalParameterFunction;
        std::shared_ptr< GravityFieldModel > gravityFieldCentralBody = bodyExertingAcceleration->getGravityFieldModel( );
        if( gravityFieldCentralBody == nullptr )
        {
            throw std::runtime_error( "Error " + nameOfBodyExertingAcceleration + " does not have a gravity field " +
                                      "when making SEP violation acceleration on" + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            centralBodyGravitationalParameterFunction =
                    std::bind( &GravityFieldModel::getGravitationalParameter, bodyExertingAcceleration->getGravityFieldModel( ) );
        }

        // get SEP corrected position of central body (see functions above)
        std::function< Eigen::Vector3d( ) > sepPositionCorrectionFunction =
                std::bind( getSEPCorrectedPosition,
                                bodyExertingAcceleration,
                                nameOfBodyExertingAcceleration,
                                bodyMap,
                                sepViolationAccelerationSettings->bodyNames_,
                                nordtvedtParameterFunction);

        // get Nordtvedt partial (see functions above)
        std::function< Eigen::Vector3d( ) > nordtvedtPartialFunction =
                std::bind( getNordtvedtPartial,
                                bodyExertingAcceleration,
                                bodyUndergoingAcceleration,
                                nameOfBodyExertingAcceleration,
                                nameOfBodyUndergoingAcceleration,
                                bodyMap,
                                sepViolationAccelerationSettings->bodyNames_);

        accelerationModel = std::make_shared< SEPViolationAcceleration >
                ( positionFunctionOfBodyUndergoingAcceleration,
                  positionFunctionOfBodyExertingAcceleration,
                  sepPositionCorrectionFunction,
                  centralBodyGravitationalParameterFunction,
                  nordtvedtPartialFunction,
                  useNordtvedtConstraintFunction,
                  nordtvedtParameterFunction
                  );

//        std::cout<<nordtvedtParameterFunction()<<" "<<sepPositionCorrectionFunction().transpose()<<std::endl;


    }
    return accelerationModel;
}

//! Function to create an orbiter relativistic correction acceleration model
std::shared_ptr< relativity::RelativisticAccelerationCorrection > createRelativisticCorrectionAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const NamedBodyMap& bodyMap )
{
    using namespace relativity;

    // Declare pointer to return object
    std::shared_ptr< RelativisticAccelerationCorrection > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< RelativisticAccelerationCorrectionSettings > relativisticAccelerationSettings =
            std::dynamic_pointer_cast< RelativisticAccelerationCorrectionSettings >(
                accelerationSettings );
    if( relativisticAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error, expected relativistic acceleration settings when making acceleration model on " +
                                  nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration );
    }
    else
    {

        // Retrieve function pointers for properties of bodies exerting/undergoing acceleration.
        std::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingAcceleration =
                std::bind( &Body::getState, bodyExertingAcceleration );
        std::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingAcceleration =
                std::bind( &Body::getState, bodyUndergoingAcceleration );

        std::function< double( ) > centralBodyGravitationalParameterFunction;
        std::shared_ptr< GravityFieldModel > gravityField = bodyExertingAcceleration->getGravityFieldModel( );
        if( gravityField == nullptr )
        {
            throw std::runtime_error( "Error " + nameOfBodyExertingAcceleration + " does not have a gravity field " +
                                      "when making relativistic acceleration on" + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            centralBodyGravitationalParameterFunction =
                    std::bind( &GravityFieldModel::getGravitationalParameter, bodyExertingAcceleration->getGravityFieldModel( ) );
        }

        std::function< double( ) > acceleratedBodyGravitationalParameterFunction;
        std::shared_ptr< GravityFieldModel > gravityField2 = bodyUndergoingAcceleration->getGravityFieldModel( );
        if( gravityField2 == nullptr )
        {
            std::cout<<"Warning: " + nameOfBodyUndergoingAcceleration + " does not have a gravity field " +
                                      "when making relativistic acceleration caused by" + nameOfBodyExertingAcceleration <<std::endl;
            std::cout<<"Using gravitational parameter = 0"<<std::endl;
            acceleratedBodyGravitationalParameterFunction = [ = ]( ){ return 0.0; };
        }
        else
        {
            acceleratedBodyGravitationalParameterFunction =
                    std::bind( &GravityFieldModel::getGravitationalParameter, bodyUndergoingAcceleration->getGravityFieldModel( ) );
        }



        std::function< double( ) > ppnGammaFunction = std::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet );
        std::function< double( ) > ppnBetaFunction = std::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet );
        std::function< double( ) > ppnAlpha1Function = std::bind( &PPNParameterSet::getParameterAlpha1, ppnParameterSet );
        std::function< double( ) > ppnAlpha2Function = std::bind( &PPNParameterSet::getParameterAlpha2, ppnParameterSet );

        // Create acceleration model if only schwarzschild term is to be used.
        if( relativisticAccelerationSettings->calculateLenseThirringCorrection_ == false &&
                relativisticAccelerationSettings->calculateDeSitterCorrection_ == false )
        {
            // Create acceleration model.
            accelerationModel = std::make_shared< RelativisticAccelerationCorrection >
                    ( stateFunctionOfBodyUndergoingAcceleration,
                      stateFunctionOfBodyExertingAcceleration,
                      centralBodyGravitationalParameterFunction,
                      acceleratedBodyGravitationalParameterFunction,
                      ppnGammaFunction, ppnBetaFunction,
                      ppnAlpha1Function, ppnAlpha2Function);

        }
        else
        {

            // Retrieve parameters of primary body if de Sitter term is to be used.
            std::function< Eigen::Vector6d( ) > stateFunctionOfPrimaryBody;
            std::function< double( ) > primaryBodyGravitationalParameterFunction;
            if( relativisticAccelerationSettings->calculateDeSitterCorrection_ == true )
            {
                if(  bodyMap.count( relativisticAccelerationSettings->primaryBody_ ) == 0 )
                {
                    throw std::runtime_error( "Error, no primary body " + relativisticAccelerationSettings->primaryBody_ +
                                              " found when making de Sitter acceleration correction" );
                }
                stateFunctionOfPrimaryBody =
                        std::bind( &Body::getState, bodyMap.at( relativisticAccelerationSettings->primaryBody_ ) );

                if(  bodyMap.at( relativisticAccelerationSettings->primaryBody_ )->getGravityFieldModel( ) == nullptr )
                {
                    throw std::runtime_error( "Error, primary body " + relativisticAccelerationSettings->primaryBody_ +
                                              " has no gravity field when making de Sitter acceleration correction" );
                }

                primaryBodyGravitationalParameterFunction =
                        std::bind( &GravityFieldModel::getGravitationalParameter,
                                   bodyMap.at( relativisticAccelerationSettings->primaryBody_ )->getGravityFieldModel( ) );


            }

            // Retrieve angular momentum vector if Lense-Thirring
            if( relativisticAccelerationSettings->calculateLenseThirringCorrection_ == true  ){

                // check if rotation model is set
                std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris = bodyExertingAcceleration->getRotationalEphemeris( );
                if( rotationalEphemeris  == nullptr )
                {
                    throw std::runtime_error( "Error " + nameOfBodyExertingAcceleration + " needs to have a simple rotational ephemeris set to calculate Lense-Thirring acceleration ");
                }

                if( rotationalEphemeris->getAngularMomentum() == 0.0){
                    throw std::runtime_error( "Error " + nameOfBodyExertingAcceleration + " has an angular momentum of 0, asked Lense-Thirring acceleration does not exist ");
                }


                // get constant angular momentum along z-axis of central body
                std::function< double( ) > angularMomentumFunction;
//                angularMomentumFunction = [ = ]( ){ return
//                            relativisticAccelerationSettings->centralBodyAngularMomentum_;};
                angularMomentumFunction =
                            std::bind( &ephemerides::RotationalEphemeris::getAngularMomentum,
                                       bodyExertingAcceleration->getRotationalEphemeris( ));


                // get transformation quaternoind from body local frame to global frame
                std::function< Eigen::Quaterniond( const double ) > quaternoidFromCentralBodyToGlobalFrameFunction;
                quaternoidFromCentralBodyToGlobalFrameFunction =
                    std::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame,
                               bodyExertingAcceleration->getRotationalEphemeris( ),
                               std::placeholders::_1 );

                // function to retrieve central body name
                std::function< std::string( ) > nameOfBodyExertingAccelerationFunction;
                nameOfBodyExertingAccelerationFunction = [ = ]( ){ return
                            nameOfBodyExertingAcceleration; };

                if( relativisticAccelerationSettings->calculateDeSitterCorrection_ == true )
                {
                    // Create acceleration model with Lense-Thirring and de Sitter terms.
                    accelerationModel = std::make_shared< RelativisticAccelerationCorrection >
                            ( stateFunctionOfBodyUndergoingAcceleration,
                              stateFunctionOfBodyExertingAcceleration,
                              stateFunctionOfPrimaryBody,
                              centralBodyGravitationalParameterFunction,
                              acceleratedBodyGravitationalParameterFunction,
                              primaryBodyGravitationalParameterFunction,
                              relativisticAccelerationSettings->primaryBody_,
                              angularMomentumFunction,
                              quaternoidFromCentralBodyToGlobalFrameFunction,
                              nameOfBodyExertingAccelerationFunction,
                              ppnGammaFunction, ppnBetaFunction,
                              ppnAlpha1Function, ppnAlpha2Function,
                              relativisticAccelerationSettings->calculateSchwarzschildCorrection_ );
                }
                else
                {

                    // Create acceleration model with Lense-Thirring and term.
                    accelerationModel = std::make_shared< RelativisticAccelerationCorrection >
                            ( stateFunctionOfBodyUndergoingAcceleration,
                              stateFunctionOfBodyExertingAcceleration,
                              centralBodyGravitationalParameterFunction,
                              acceleratedBodyGravitationalParameterFunction,
                              angularMomentumFunction,
                              quaternoidFromCentralBodyToGlobalFrameFunction,
                              nameOfBodyExertingAccelerationFunction,
                              ppnGammaFunction, ppnBetaFunction,
                              ppnAlpha1Function, ppnAlpha2Function,
                              relativisticAccelerationSettings->calculateSchwarzschildCorrection_ );
                }
            }
        }
    }
    return accelerationModel;
}


//! Function to create empirical acceleration model.
std::shared_ptr< EmpiricalAcceleration > createEmpiricalAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Declare pointer to return object
    std::shared_ptr< EmpiricalAcceleration > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< EmpiricalAccelerationSettings > empiricalSettings =
            std::dynamic_pointer_cast< EmpiricalAccelerationSettings >(
                accelerationSettings );
    if( empiricalSettings == nullptr )
    {
        throw std::runtime_error( "Error, expected empirical acceleration settings when making acceleration model on " +
                                  nameOfBodyUndergoingAcceleration + " due to " + nameOfBodyExertingAcceleration );
    }
    else
    {
        // Get pointer to gravity field of central body (for determining keplerian elememts)
        std::shared_ptr< GravityFieldModel > gravityField = bodyExertingAcceleration->getGravityFieldModel( );

        if( gravityField == nullptr )
        {
            throw std::runtime_error( "Error " + nameOfBodyExertingAcceleration + " does not have a gravity field " +
                                      "when making empirical acceleration on" + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            // Create acceleration model.
            accelerationModel = std::make_shared< EmpiricalAcceleration >(
                        empiricalSettings->constantAcceleration_,
                        empiricalSettings->sineAcceleration_,
                        empiricalSettings->cosineAcceleration_,
                        std::bind( &Body::getState, bodyUndergoingAcceleration ),
                        std::bind( &GravityFieldModel::getGravitationalParameter, gravityField ),
                        std::bind( &Body::getState, bodyExertingAcceleration ) );
        }
    }

    return accelerationModel;
}

//! Function to create a thrust acceleration model.
std::shared_ptr< propulsion::ThrustAcceleration >
createThrustAcceleratioModel(
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyUndergoingThrust )
{
    // Check input consistency
    std::shared_ptr< ThrustAccelerationSettings > thrustAccelerationSettings =
            std::dynamic_pointer_cast< ThrustAccelerationSettings >( accelerationSettings );
    if( thrustAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating thrust acceleration, input is inconsistent" );
    }

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > magnitudeUpdateSettings;
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > directionUpdateSettings;



    // Check if user-supplied interpolator for full thrust ius present.
    if( thrustAccelerationSettings->interpolatorInterface_ != nullptr )
    {
        // Check input consisten
        if( thrustAccelerationSettings->thrustFrame_ == unspecified_thurst_frame )
        {
            throw std::runtime_error( "Error when creating thrust acceleration, input frame is inconsistent with interface" );
        }
        else if( thrustAccelerationSettings->thrustFrame_ != inertial_thurst_frame )
        {
            // Create rotation function from thrust-frame to propagation frame.
            if( thrustAccelerationSettings->thrustFrame_ == lvlh_thrust_frame )
            {
                std::function< Eigen::Vector6d( ) > vehicleStateFunction =
                        std::bind( &Body::getState, bodyMap.at( nameOfBodyUndergoingThrust ) );
                std::function< Eigen::Vector6d( ) > centralBodyStateFunction;

                if( ephemerides::isFrameInertial( thrustAccelerationSettings->centralBody_ ) )
                {
                    centralBodyStateFunction =  [ ]( ){ return Eigen::Vector6d::Zero( ); };
                }
                else
                {
                    if( bodyMap.count( thrustAccelerationSettings->centralBody_ ) == 0 )
                    {
                        throw std::runtime_error( "Error when creating thrust acceleration, input central body not found" );
                    }
                    centralBodyStateFunction =
                            std::bind( &Body::getState, bodyMap.at( thrustAccelerationSettings->centralBody_ ) );
                }
                thrustAccelerationSettings->interpolatorInterface_->resetRotationFunction(
                            std::bind( &reference_frames::getVelocityBasedLvlhToInertialRotationFromFunctions,
                                       vehicleStateFunction, centralBodyStateFunction, true ) );
            }
            else
            {
                throw std::runtime_error( "Error when creating thrust acceleration, input frame not recognized" );
            }
        }
    }

    // Create thrust direction model.
    std::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > thrustDirectionGuidance = createThrustGuidanceModel(
                thrustAccelerationSettings->thrustDirectionGuidanceSettings_, bodyMap, nameOfBodyUndergoingThrust,
                getBodyFixedThrustDirection( thrustAccelerationSettings->thrustMagnitudeSettings_, bodyMap,
                                             nameOfBodyUndergoingThrust ), magnitudeUpdateSettings );

    // Create thrust magnitude model
    std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitude = createThrustMagnitudeWrapper(
                thrustAccelerationSettings->thrustMagnitudeSettings_, bodyMap, nameOfBodyUndergoingThrust,
                directionUpdateSettings );

    // Add required updates of environemt models.
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > totalUpdateSettings;
    propagators::addEnvironmentUpdates( totalUpdateSettings, magnitudeUpdateSettings );
    propagators::addEnvironmentUpdates( totalUpdateSettings, directionUpdateSettings );

    // Set DependentOrientationCalculator for body if required.
    if( !( thrustAccelerationSettings->thrustDirectionGuidanceSettings_->thrustDirectionType_ ==
           thrust_direction_from_existing_body_orientation ) )
    {
        bodyMap.at( nameOfBodyUndergoingThrust )->setDependentOrientationCalculator( thrustDirectionGuidance );
    }

    // Create and return thrust acceleration object.
    std::function< void( const double ) > updateFunction =
            std::bind( &updateThrustMagnitudeAndDirection, thrustMagnitude, thrustDirectionGuidance, std::placeholders::_1 );
    std::function< void( const double ) > timeResetFunction =
            std::bind( &resetThrustMagnitudeAndDirectionTime, thrustMagnitude, thrustDirectionGuidance, std::placeholders::_1 );
    return std::make_shared< propulsion::ThrustAcceleration >(
                std::bind( &propulsion::ThrustMagnitudeWrapper::getCurrentThrustMagnitude, thrustMagnitude ),
                std::bind( &propulsion::BodyFixedForceDirectionGuidance ::getCurrentForceDirectionInPropagationFrame, thrustDirectionGuidance ),
                std::bind( &Body::getBodyMass, bodyMap.at( nameOfBodyUndergoingThrust ) ),
                std::bind( &propulsion::ThrustMagnitudeWrapper::getCurrentMassRate, thrustMagnitude ),
                thrustAccelerationSettings->thrustMagnitudeSettings_->thrustOriginId_,
                updateFunction, timeResetFunction, totalUpdateSettings );
}

//! Function to create a direct tical acceleration model, according to approach of Lainey et al. (2007, 2009, ...)
std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > createDirectTidalDissipationAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Check input consistency
    std::shared_ptr< DirectTidalDissipationAccelerationSettings > tidalAccelerationSettings =
            std::dynamic_pointer_cast< DirectTidalDissipationAccelerationSettings >( accelerationSettings );
    if( tidalAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, input is inconsistent" );
    }

    std::function< double( ) > gravitationalParaterFunctionOfBodyExertingTide;
    std::function< double( ) > gravitationalParaterFunctionOfBodyUndergoingTide;

    if( tidalAccelerationSettings->useTideRaisedOnPlanet_ )
    {
        if( bodyUndergoingAcceleration->getGravityFieldModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, satellite " +
                                      nameOfBodyUndergoingAcceleration + " has no gravity field" );
        }
        else
        {
            gravitationalParaterFunctionOfBodyUndergoingTide = std::bind(
                        &GravityFieldModel::getGravitationalParameter, bodyUndergoingAcceleration->getGravityFieldModel( ) );
        }
    }
    else
    {
        if( bodyExertingAcceleration->getGravityFieldModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, satellite " +
                                      nameOfBodyExertingAcceleration + " has no gravity field" );
        }
        else
        {
            gravitationalParaterFunctionOfBodyExertingTide = std::bind(
                        &GravityFieldModel::getGravitationalParameter, bodyExertingAcceleration->getGravityFieldModel( ) );
        }


        if( bodyUndergoingAcceleration->getGravityFieldModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, satellite " +
                                      nameOfBodyUndergoingAcceleration + " has no gravity field" );
        }
        else
        {
            gravitationalParaterFunctionOfBodyUndergoingTide = std::bind(
                        &GravityFieldModel::getGravitationalParameter, bodyUndergoingAcceleration->getGravityFieldModel( ) );
        }
    }

    double referenceRadius = TUDAT_NAN;
    if( tidalAccelerationSettings->useTideRaisedOnPlanet_ )
    {
        if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) ) == nullptr )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, planet " +
                                      nameOfBodyExertingAcceleration + " has no s.h. gravity field" );
        }
        else
        {
            referenceRadius = std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                        bodyExertingAcceleration->getGravityFieldModel( ) )->getReferenceRadius( );
        }
    }
    else
    {
        if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                    bodyUndergoingAcceleration->getGravityFieldModel( ) ) == nullptr )
        {
            throw std::runtime_error( "Error when creating direct tidal dissipation acceleration, planet " +
                                      nameOfBodyUndergoingAcceleration + " has no s.h. gravity field" );
        }
        else
        {
            referenceRadius = std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                        bodyUndergoingAcceleration->getGravityFieldModel( ) )->getReferenceRadius( );
        }
    }
    

    if( tidalAccelerationSettings->useTideRaisedOnPlanet_ )
    {
        std::function< Eigen::Vector3d( ) > planetAngularVelocityVectorFunction =
                std::bind( &Body::getCurrentAngularVelocityVectorInGlobalFrame, bodyExertingAcceleration );


        return std::make_shared< DirectTidalDissipationAcceleration >(
                    std::bind( &Body::getState, bodyUndergoingAcceleration ),
                    std::bind( &Body::getState, bodyExertingAcceleration ),
                    gravitationalParaterFunctionOfBodyUndergoingTide,
                    planetAngularVelocityVectorFunction,
                    tidalAccelerationSettings->k2LoveNumber_,
                    tidalAccelerationSettings->timeLag_,
                    referenceRadius,
                    tidalAccelerationSettings->includeDirectRadialComponent_);
    }
    else
    {
        return std::make_shared< DirectTidalDissipationAcceleration >(
                    std::bind( &Body::getState, bodyUndergoingAcceleration ),
                    std::bind( &Body::getState, bodyExertingAcceleration ),
                    gravitationalParaterFunctionOfBodyExertingTide,
                    gravitationalParaterFunctionOfBodyUndergoingTide,
                    tidalAccelerationSettings->k2LoveNumber_,
                    tidalAccelerationSettings->timeLag_,
                    referenceRadius,
                    tidalAccelerationSettings->includeDirectRadialComponent_);
    }
}

//! Function to create a momentum wheel desaturation acceleration model.
std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > createMomentumWheelDesaturationAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  std::shared_ptr< AccelerationSettings > accelerationSettings )
{
    // Check input consistency
    std::shared_ptr< MomentumWheelDesaturationAccelerationSettings > desaturationAccelerationSettings =
            std::dynamic_pointer_cast< MomentumWheelDesaturationAccelerationSettings >( accelerationSettings );
    if( desaturationAccelerationSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating momentum wheel desaturation acceleration, input is inconsistent" );
    }

    if( nameOfBodyUndergoingAcceleration != nameOfBodyExertingAcceleration )
    {
        throw std::runtime_error( "Error when creating momentum wheel desaturation acceleration, exerting and undergoing bodies are not the same" );
    }

    // Return desaturation acceleration model.
    return std::make_shared< propulsion::MomentumWheelDesaturationThrustAcceleration >(
                desaturationAccelerationSettings->thrustMidTimes_,
                desaturationAccelerationSettings->deltaVValues_,
                desaturationAccelerationSettings->totalManeuverTime_,
                desaturationAccelerationSettings->maneuverRiseTime_ );
}

//! Function to create acceleration model object.
std::shared_ptr< AccelerationModel< Eigen::Vector3d > > createAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfCentralBody,
        const NamedBodyMap& bodyMap )
{
    // Declare pointer to return object.
    std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelPointer;

    // Switch to call correct acceleration model type factory function.
    switch( accelerationSettings->accelerationType_ )
    {
    case central_gravity:
        accelerationModelPointer = createGravitationalAccelerationModel(
                    bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationSettings,
                    nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                    centralBody, nameOfCentralBody );
        break;
    case spherical_harmonic_gravity:
        accelerationModelPointer = createGravitationalAccelerationModel(
                    bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationSettings,
                    nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                    centralBody, nameOfCentralBody );
        break;
    case mutual_spherical_harmonic_gravity:
        accelerationModelPointer = createGravitationalAccelerationModel(
                    bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationSettings,
                    nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
                    centralBody, nameOfCentralBody );
        break;
    case aerodynamic:
        accelerationModelPointer = createAerodynamicAcceleratioModel(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration );
        break;
    case cannon_ball_radiation_pressure:
        accelerationModelPointer = createCannonballRadiationPressureAcceleratioModel(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration );
        break;
    case panelled_radiation_pressure_acceleration:
        accelerationModelPointer = createPanelledRadiationPressureAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration );
        break;
    case thrust_acceleration:
        accelerationModelPointer = createThrustAcceleratioModel(
                    accelerationSettings, bodyMap,
                    nameOfBodyUndergoingAcceleration );
        break;
    case relativistic_correction_acceleration:
        accelerationModelPointer = createRelativisticCorrectionAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings, bodyMap );
        break;
    case time_varying_gravitational_parameter_acceleration:
        accelerationModelPointer = createTimeVaryingGravitationalParameterAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings);
        break;
    case sep_violation_acceleration:
        accelerationModelPointer = createSEPViolationAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings,
                    bodyMap);
        break;
    case empirical_acceleration:
        accelerationModelPointer = createEmpiricalAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    case direct_tidal_dissipation_in_central_body_acceleration:
        accelerationModelPointer = createDirectTidalDissipationAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    case direct_tidal_dissipation_in_orbiting_body_acceleration:
        accelerationModelPointer = createDirectTidalDissipationAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    case momentum_wheel_desaturation_acceleration:
        accelerationModelPointer = createMomentumWheelDesaturationAcceleration(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration,
                    accelerationSettings );
        break;
    case solar_sail_acceleration:
        accelerationModelPointer = createSolarSailAccelerationModel(
                    bodyUndergoingAcceleration,
                    bodyExertingAcceleration,
                    centralBody,
                    nameOfBodyUndergoingAcceleration,
                    nameOfBodyExertingAcceleration );
        break;
    default:
        throw std::runtime_error(
                    std::string( "Error, acceleration model ") +
                    std::to_string( accelerationSettings->accelerationType_ ) +
                    " not recognized when making acceleration model of" +
                    nameOfBodyExertingAcceleration + " on " +
                    nameOfBodyUndergoingAcceleration );
        break;
    }
    return accelerationModelPointer;
}

//! Function to put SelectedAccelerationMap in correct order, to ensure correct model creation
SelectedAccelerationList orderSelectedAccelerationMap( const SelectedAccelerationMap& selectedAccelerationsPerBody )
{
    // Declare map of acceleration models acting on current body.
    SelectedAccelerationList orderedAccelerationsPerBody;

    // Iterate over all bodies which are undergoing acceleration
    for( SelectedAccelerationMap::const_iterator bodyIterator =
         selectedAccelerationsPerBody.begin( ); bodyIterator != selectedAccelerationsPerBody.end( );
         bodyIterator++ )
    {
        // Retrieve name of body undergoing acceleration.
        std::string bodyUndergoingAcceleration = bodyIterator->first;

        // Retrieve list of required acceleration model types and bodies exerting accelerationd on
        // current body.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > >
                accelerationsForBody = bodyIterator->second;

        // Retrieve indices of all acceleration anf thrust models.
        std::vector< int > aerodynamicAccelerationIndices;
        std::vector< int > thrustAccelerationIndices;

        std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > >
                currentBodyAccelerations;
        int counter = 0;
        // Iterate over all bodies exerting an acceleration
        for( std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > >::
             iterator body2Iterator = accelerationsForBody.begin( );
             body2Iterator != accelerationsForBody.end( ); body2Iterator++ )
        {
            // Retrieve name of body exerting acceleration.
            std::string bodyExertingAcceleration = body2Iterator->first;
            std::vector< std::shared_ptr< AccelerationSettings > > accelerationList = body2Iterator->second;
            for( unsigned int i = 0; i < accelerationList.size( ); i++ )
            {
                if( accelerationList.at( i )->accelerationType_ == basic_astrodynamics::thrust_acceleration )
                {
                    thrustAccelerationIndices.push_back( counter );
                }
                else if( accelerationList.at( i )->accelerationType_ == basic_astrodynamics::aerodynamic )
                {
                    aerodynamicAccelerationIndices.push_back( counter );
                }
                std::pair< std::string, std::shared_ptr< AccelerationSettings > >  currentAccelerationPair =
                        std::make_pair( bodyExertingAcceleration, accelerationList.at( i ) );
                currentBodyAccelerations.push_back( currentAccelerationPair );
                counter++;
            }
        }

        if( thrustAccelerationIndices.size( ) > 0 && aerodynamicAccelerationIndices.size( ) > 0 )
        {
            std::vector< int > indexList;
            for( unsigned int i = 0; i < aerodynamicAccelerationIndices.size( ); i++ )
            {
                indexList.push_back( aerodynamicAccelerationIndices.at( i ) );
            }
            for( unsigned int i = 0; i < thrustAccelerationIndices.size( ); i++ )
            {
                indexList.push_back( thrustAccelerationIndices.at( i ) );
            }

            std::vector< int > unorderedIndexList = indexList;
            std::sort( indexList.begin( ), indexList.end( ) );
            if( !( indexList == unorderedIndexList ) )
            {
                std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > >
                        orderedAccelerationSettings = currentBodyAccelerations;

                int indexCounter = 0;
                for( unsigned int i = 0; i < aerodynamicAccelerationIndices.size( ); i++ )
                {
                    orderedAccelerationSettings[ indexList.at( indexCounter ) ]
                            = currentBodyAccelerations[ aerodynamicAccelerationIndices.at( i ) ];
                    indexCounter++;
                }

                for( unsigned int i = 0; i < thrustAccelerationIndices.size( ); i++ )
                {
                    orderedAccelerationSettings[ indexList.at( indexCounter ) ]
                            = currentBodyAccelerations[ thrustAccelerationIndices.at( i ) ];
                    indexCounter++;
                }

                currentBodyAccelerations = orderedAccelerationSettings;
            }
        }

        orderedAccelerationsPerBody[ bodyUndergoingAcceleration ] = currentBodyAccelerations;
    }

    return orderedAccelerationsPerBody;
}


//! Function to create a set of acceleration models from a map of bodies and acceleration model types.
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::map< std::string, std::string >& centralBodies )
{
    // Declare return map.
    basic_astrodynamics::AccelerationMap accelerationModelMap;

    // Put selectedAccelerationPerBody in correct order
    SelectedAccelerationList orderedAccelerationPerBody =
            orderSelectedAccelerationMap( selectedAccelerationPerBody );

    // Iterate over all bodies which are undergoing acceleration
    for( SelectedAccelerationList::const_iterator bodyIterator =
         orderedAccelerationPerBody.begin( ); bodyIterator != orderedAccelerationPerBody.end( );
         bodyIterator++ )
    {
        std::shared_ptr< Body > currentCentralBody;

        // Retrieve name of body undergoing acceleration.
        std::string bodyUndergoingAcceleration = bodyIterator->first;

        // Retrieve name of current central body.
        std::string currentCentralBodyName = centralBodies.at( bodyUndergoingAcceleration );

        if( !ephemerides::isFrameInertial( currentCentralBodyName ) )
        {
            if( bodyMap.count( currentCentralBodyName ) == 0 )
            {
                throw std::runtime_error(
                            std::string( "Error, could not find non-inertial central body ") +
                            currentCentralBodyName + " of " + bodyUndergoingAcceleration +
                            " when making acceleration model." );
            }
            else
            {
                currentCentralBody = bodyMap.at( currentCentralBodyName );
            }
        }

        // Check if body undergoing acceleration is included in bodyMap
        if( bodyMap.count( bodyUndergoingAcceleration ) ==  0 )
        {
            throw std::runtime_error(
                        std::string( "Error when making acceleration models, requested forces" ) +
                        "acting on body " + bodyUndergoingAcceleration  +
                        ", but no such body found in map of bodies" );
        }

        // Declare map of acceleration models acting on current body.
        basic_astrodynamics::SingleBodyAccelerationMap mapOfAccelerationsForBody;

        // Retrieve list of required acceleration model types and bodies exerting accelerationd on
        // current body.
        std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > >
                accelerationsForBody = bodyIterator->second;

        std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > > thrustAccelerationSettings;

        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > currentAcceleration;
        // Iterate over all bodies exerting an acceleration
        for( unsigned int i = 0; i < accelerationsForBody.size( ); i++ )
        {
            // Retrieve name of body exerting acceleration.
            std::string bodyExertingAcceleration = accelerationsForBody.at( i ).first;

            // Check if body exerting acceleration is included in bodyMap
            if( bodyMap.count( bodyExertingAcceleration ) ==  0 )
            {
                throw std::runtime_error(
                            std::string( "Error when making acceleration models, requested forces ")
                            + "acting on body " + bodyUndergoingAcceleration  + " due to body " +
                            bodyExertingAcceleration +
                            ", but no such body found in map of bodies" );
            }

            if( !( accelerationsForBody.at( i ).second->accelerationType_ == basic_astrodynamics::thrust_acceleration ) )
            {
                currentAcceleration = createAccelerationModel( bodyMap.at( bodyUndergoingAcceleration ),
                                                               bodyMap.at( bodyExertingAcceleration ),
                                                               accelerationsForBody.at( i ).second,
                                                               bodyUndergoingAcceleration,
                                                               bodyExertingAcceleration,
                                                               currentCentralBody,
                                                               currentCentralBodyName,
                                                               bodyMap );


                // Create acceleration model.
                mapOfAccelerationsForBody[ bodyExertingAcceleration ].push_back(
                            currentAcceleration );
            }
            else
            {
                thrustAccelerationSettings.push_back( accelerationsForBody.at( i ) );
            }

        }

        for( unsigned int i = 0; i < thrustAccelerationSettings.size( ); i++ )
        {
            currentAcceleration = createAccelerationModel( bodyMap.at( bodyUndergoingAcceleration ),
                                                           bodyMap.at( thrustAccelerationSettings.at( i ).first ),
                                                           thrustAccelerationSettings.at( i ).second,
                                                           bodyUndergoingAcceleration,
                                                           thrustAccelerationSettings.at( i ).first,
                                                           currentCentralBody,
                                                           currentCentralBodyName,
                                                           bodyMap );


            // Create acceleration model.
            mapOfAccelerationsForBody[ thrustAccelerationSettings.at( i ).first  ].push_back(
                        currentAcceleration );
        }


        // Put acceleration models on current body in return map.
        accelerationModelMap[ bodyUndergoingAcceleration ] = mapOfAccelerationsForBody;
    }

    return accelerationModelMap;
}

//! Function to create acceleration models from a map of bodies and acceleration model types.
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::vector< std::string >& propagatedBodies,
        const std::vector< std::string >& centralBodies )
{
    if( centralBodies.size( ) != propagatedBodies.size( ) )
    {
        throw std::runtime_error( "Error, number of propagated bodies must equal number of central bodies" );
    }

    std::map< std::string, std::string > centralBodyMap;
    for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
    {
        centralBodyMap[ propagatedBodies.at( i ) ] = centralBodies.at( i );
    }

    return createAccelerationModelsMap( bodyMap, selectedAccelerationPerBody, centralBodyMap );
}

} // namespace simulation_setup

} // namespace tudat
