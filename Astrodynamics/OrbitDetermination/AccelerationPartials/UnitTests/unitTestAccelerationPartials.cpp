/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantDragCoefficient.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/empiricalAccelerationCoefficients.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/gravitationalParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/radiationPressureCoefficient.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/ppnParameters.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/directTidalTimeLag.h"
#include "Tudat/Astrodynamics/Relativity/metric.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/numericalAccelerationPartial.h"
#include "Tudat/Astrodynamics/Relativity/relativisticAccelerationCorrection.h"
#include "Tudat/SimulationSetup/EstimationSetup/createAccelerationPartials.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"

#include "tudatApplications/thesis/MyApplications/sepViolationAcceleration.h"
#include "tudatApplications/thesis/MyApplications/timeVaryingGravitationalParameterAcceleration.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::relativity;
using namespace tudat::gravitation;
using namespace tudat::aerodynamics;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::unit_conversions;
using namespace tudat::orbit_determination;
using namespace tudat::acceleration_partials;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::electro_magnetism;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_acceleration_partials )

BOOST_AUTO_TEST_CASE( testCentralGravityPartials )
{
    // Create empty bodies, earth and sun.
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    std::shared_ptr< Body > sun = std::make_shared< Body >( );

    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = earth;
    bodyMap[ "Sun" ] = sun;

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "Sun", "J2000", "NONE", 1.0E6 ) );
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "Sun", "J2000", "NONE", 1.0E6 ) );

    // Get sun gravitational parameter and set gravity field model.
    double sunsGravitationalParameter = getBodyGravitationalParameter( "Sun" );
    std::shared_ptr< GravityFieldModel > sunGravityFieldModel =
            std::make_shared< GravityFieldModel >( sunsGravitationalParameter );
    sun->setGravityFieldModel( sunGravityFieldModel );
    double earthGravitationalParameter = getBodyGravitationalParameter( "Earth" );
    std::shared_ptr< GravityFieldModel > earthGravityFieldModel =
            std::make_shared< GravityFieldModel >( earthGravitationalParameter );
    earth->setGravityFieldModel( earthGravityFieldModel );

    // Create acceleration due to sun on earth.
    std::shared_ptr< CentralGravitationalAccelerationModel3d > gravitationalAcceleration =\
            createCentralGravityAcceleratioModel( earth, sun, "Earth", "Sun", 1 );

    // Create central gravity partial.
    std::shared_ptr< AccelerationPartial > centralGravitationPartial =
            createAnalyticalAccelerationPartial( gravitationalAcceleration, std::make_pair( "Earth", earth ),
                                                 std::make_pair( "Sun", sun ), bodyMap );

    // Create gravitational parameter object.
    std::shared_ptr< EstimatableParameter< double > > sunGravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( sunGravityFieldModel, "Sun" );
    std::shared_ptr< EstimatableParameter< double > > earthGravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( earthGravityFieldModel, "Earth" );

    // Calculate analytical partials.
    centralGravitationPartial->update( 0.0 );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    centralGravitationPartial->wrtPositionOfAcceleratedBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    centralGravitationPartial->wrtVelocityOfAcceleratedBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
    centralGravitationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    centralGravitationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::Vector3d partialWrtSunGravitationalParameter = centralGravitationPartial->wrtParameter(
                sunGravitationalParameterParameter );
    Eigen::Vector3d partialWrtEarthGravitationalParameter = centralGravitationPartial->wrtParameter(
                earthGravitationalParameterParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10000.0, 10000.0, 10000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, earth, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            std::bind( &Body::setState, sun, std::placeholders::_1 );
    std::function< Eigen::Vector6d ( ) > earthStateGetFunction =
            std::bind( &Body::getState, earth );
    std::function< Eigen::Vector6d ( ) > sunStateGetFunction =
            std::bind( &Body::getState, sun );

    // Calculate numerical partials.
    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), positionPerturbation, 0 );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), velocityPerturbation, 3 );
    testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, gravitationalAcceleration, sun->getState( ), positionPerturbation, 0 );
    testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, gravitationalAcceleration, sun->getState( ), velocityPerturbation, 3 );
    Eigen::Vector3d testPartialWrtSunGravitationalParameter = calculateAccelerationWrtParameterPartials(
                sunGravitationalParameterParameter, gravitationalAcceleration, 1.0E12 );
    Eigen::Vector3d testPartialWrtEarthGravitationalParameter = calculateAccelerationWrtParameterPartials(
                earthGravitationalParameterParameter, gravitationalAcceleration, 1.0E12 );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                       partialWrtEarthPosition, 1.0E-8 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                       partialWrtEarthVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                       partialWrtSunPosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                       partialWrtSunVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunGravitationalParameter,
                                       partialWrtSunGravitationalParameter, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtEarthGravitationalParameter,
                                       testPartialWrtEarthGravitationalParameter, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtEarthGravitationalParameter,
                                       partialWrtSunGravitationalParameter, std::numeric_limits< double >::epsilon(  ) );
}

BOOST_AUTO_TEST_CASE( testRadiationPressureAccelerationPartials )
{
    // Create empty bodies, earth and sun.
    std::shared_ptr< Body > vehicle = std::make_shared< Body >( );
    vehicle->setConstantBodyMass( 400.0 );
    std::shared_ptr< Body > sun = std::make_shared< Body >( );

    NamedBodyMap bodyMap;
    bodyMap[ "Vehicle" ] = vehicle;
    bodyMap[ "Sun" ] = sun;

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "SSB", "J2000", "NONE", 1.0E6 ) );
    vehicle->setState( getBodyCartesianStateAtEpoch(  "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            std::bind( &Body::setState, sun, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            std::bind( &Body::setState, vehicle, std::placeholders::_1 );
    std::function< Eigen::Vector6d( ) > sunStateGetFunction =
            std::bind( &Body::getState, sun );
    std::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
            std::bind( &Body::getState, vehicle );

    // Create radiation pressure properties of vehicle
    std::shared_ptr< RadiationPressureInterface > radiationPressureInterface =
            createRadiationPressureInterface( std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                                                  "Sun", mathematical_constants::PI * 0.3 * 0.3, 1.2 ), "Vehicle", bodyMap );
    radiationPressureInterface->updateInterface( 0.0 );
    vehicle->setRadiationPressureInterface( "Sun", radiationPressureInterface );

    // Create acceleration model.
    std::shared_ptr< CannonBallRadiationPressureAcceleration > accelerationModel =
            std::make_shared< CannonBallRadiationPressureAcceleration >(
                std::bind( &Body::getPosition, sun ),
                std::bind( &Body::getPosition, vehicle ),
                std::bind( &RadiationPressureInterface::getCurrentRadiationPressure,
                           radiationPressureInterface ),
                std::bind( &RadiationPressureInterface::getRadiationPressureCoefficient, radiationPressureInterface ),
                std::bind( &RadiationPressureInterface::getArea, radiationPressureInterface ),
                std::bind( &Body::getBodyMass, vehicle ) );

    // Create partial-calculating object.
    std::shared_ptr< AccelerationPartial > accelerationPartial =
            createAnalyticalAccelerationPartial( accelerationModel, std::make_pair( "Vehicle", vehicle ),
                                                 std::make_pair( "Sun", sun ), bodyMap );

    // Create parameter object
    std::string vehicleName = "Vehicle";
    std::shared_ptr< EstimatableParameter< double > > radiationPressureCoefficient =
            std::make_shared< RadiationPressureCoefficient >( radiationPressureInterface, vehicleName );

    std::vector< double > timeLimits;
    timeLimits.push_back( 0.0 );
    timeLimits.push_back( 3600.0 );
    timeLimits.push_back( 7200.0 );
    timeLimits.push_back( 10800.0 );

    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > arcWiseRadiationPressureCoefficient =
            std::make_shared< ArcWiseRadiationPressureCoefficient >( radiationPressureInterface, timeLimits, vehicleName );


    // Calculate analytical partials.
    double currentTime = 0.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::Vector3d partialWrtRadiationPressureCoefficient = accelerationPartial->wrtParameter(
                radiationPressureCoefficient );

    // Get arc-wise radiation pressure coefficient partials
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );
    currentTime = 1000.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise2 = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );
    currentTime = 4000.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise3 = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );
    currentTime = 7000.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise4 = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );
    currentTime = 10000.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise5 = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );
    currentTime = 12000.0;
    accelerationPartial->update( currentTime );
    Eigen::MatrixXd partialWrtRadiationPressureCoefficientArcwise6 = accelerationPartial->wrtParameter(
                arcWiseRadiationPressureCoefficient );

    // Check whether arc-wise radiation pressure partials are properly segmented
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 4; j++ )
        {
            if( j != 0 )
            {
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise2( i, j ), 0.0 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise2( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
            }

            if( j != 1 )
            {
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise3( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise4( i, j ), 0.0 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise3( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise4( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
            }

            if( j != 2 )
            {
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise5( i, j ), 0.0 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise5( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
            }

            if( j != 3 )
            {
                BOOST_CHECK_EQUAL( partialWrtRadiationPressureCoefficientArcwise6( i, j ), 0.0 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( partialWrtRadiationPressureCoefficientArcwise6( i, j ) -
                                              partialWrtRadiationPressureCoefficient( i ) ), 1.0E-24 );
            }
        }
    }

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Vector3d testPartialWrtRadiationPressureCoefficient = Eigen::Vector3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10000.0, 10000.0, 10000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Calculate numerical partials.
    std::function< void( ) > updateFunction =
            std::bind( &RadiationPressureInterface::updateInterface, radiationPressureInterface, 0.0 );
    testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, accelerationModel, sun->getState( ), positionPerturbation, 0, updateFunction );
    testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, vehicle->getState( ), positionPerturbation, 0, updateFunction );
    testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, accelerationModel, sun->getState( ),velocityPerturbation, 3, updateFunction );
    testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, vehicle->getState( ), velocityPerturbation, 3, updateFunction );
    testPartialWrtRadiationPressureCoefficient = calculateAccelerationWrtParameterPartials(
                radiationPressureCoefficient, accelerationModel, 1.0E-2, updateFunction );


    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                       partialWrtSunPosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                       partialWrtSunVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                       partialWrtVehiclePosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                       partialWrtVehicleVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtRadiationPressureCoefficient,
                                       partialWrtRadiationPressureCoefficient, 1.0E-12 );
}

BOOST_AUTO_TEST_CASE( testThirdBodyGravityPartials )
{
    // Create empty bodies, earth and sun.
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    std::shared_ptr< Body > sun = std::make_shared< Body >( );
    std::shared_ptr< Body > moon = std::make_shared< Body >( );

    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = earth;
    bodyMap[ "Sun" ] = sun;
    bodyMap[ "Moon" ] = moon;

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set current state of sun and earth.
    sun->setState( getBodyCartesianStateAtEpoch( "Sun", "Sun", "J2000", "NONE", 1.0E6 ) );
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "Sun", "J2000", "NONE", 1.0E6 ) );
    moon->setState( getBodyCartesianStateAtEpoch(  "Moon", "Sun", "J2000", "NONE", 1.0E6 ) );

    // Get sun gravitational parameter and set gravity field model.
    double sunsGravitationalParameter = getBodyGravitationalParameter( "Sun" );
    std::shared_ptr< GravityFieldModel > sunGravityFieldModel =
            std::make_shared< GravityFieldModel >( sunsGravitationalParameter );
    sun->setGravityFieldModel( sunGravityFieldModel );

    double moonsGravitationalParameter = getBodyGravitationalParameter( "Moon" );
    std::shared_ptr< GravityFieldModel > moonGravityFieldModel =
            std::make_shared< GravityFieldModel >( moonsGravitationalParameter );
    moon->setGravityFieldModel( moonGravityFieldModel );

    double earthGravitationalParameter = getBodyGravitationalParameter( "Earth" );
    std::shared_ptr< GravityFieldModel > earthGravityFieldModel =
            std::make_shared< GravityFieldModel >( earthGravitationalParameter );
    earth->setGravityFieldModel( earthGravityFieldModel );

    // Create acceleration due to moon on earth.
    std::shared_ptr< ThirdBodyCentralGravityAcceleration > gravitationalAcceleration =
            createThirdBodyCentralGravityAccelerationModel(
                moon, sun, earth, "Moon", "Sun", "Earth" );

    // Create central gravity partial.
    std::shared_ptr< AccelerationPartial > thirdBodyGravitationPartial =
            createAnalyticalAccelerationPartial( gravitationalAcceleration, std::make_pair( "Moon", moon ),
                                                 std::make_pair( "Sun", sun ), bodyMap );

    // Create gravitational parameter object.
    std::shared_ptr< EstimatableParameter< double > > gravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( sunGravityFieldModel, "Sun" );
    std::shared_ptr< EstimatableParameter< double > > moonGravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( moonGravityFieldModel, "Moon" );
    std::shared_ptr< EstimatableParameter< double > > earthGravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( earthGravityFieldModel, "Earth" );

    // Calculate analytical partials.
    thirdBodyGravitationPartial->update( 1.0E6 );
    Eigen::MatrixXd partialWrtMoonPosition = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtPositionOfAcceleratedBody( partialWrtMoonPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtMoonVelocity = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtVelocityOfAcceleratedBody( partialWrtMoonVelocity.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtPositionOfAdditionalBody( "Earth", partialWrtEarthPosition.block( 0, 0, 3, 3 )  );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    thirdBodyGravitationPartial->wrtVelocityOfAdditionalBody( "Earth", partialWrtEarthVelocity.block( 0, 0, 3, 3 )  );

    Eigen::MatrixXd partialWrtSunGravitationalParameter = thirdBodyGravitationPartial->wrtParameter(
                gravitationalParameterParameter );
    Eigen::MatrixXd partialWrtMoonGravitationalParameter = thirdBodyGravitationPartial->wrtParameter(
                moonGravitationalParameterParameter );
    Eigen::MatrixXd partialWrtEarthGravitationalParameter = thirdBodyGravitationPartial->wrtParameter(
                earthGravitationalParameterParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtMoonPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMoonVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10000.0, 10000.0, 10000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > moonStateSetFunction =
            std::bind( &Body::setState, moon, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            std::bind( &Body::setState, sun, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, earth, std::placeholders::_1 );

    // Calculate numerical partials.
    testPartialWrtMoonPosition = calculateAccelerationWrtStatePartials(
                moonStateSetFunction, gravitationalAcceleration, moon->getState( ), positionPerturbation, 0 );
    testPartialWrtMoonVelocity = calculateAccelerationWrtStatePartials(
                moonStateSetFunction, gravitationalAcceleration, moon->getState( ), velocityPerturbation, 3 );
    testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, gravitationalAcceleration, sun->getState( ), positionPerturbation, 0 );
    testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, gravitationalAcceleration, sun->getState( ), velocityPerturbation, 3 );
    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), positionPerturbation, 0 );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), velocityPerturbation, 3 );
    Eigen::Vector3d testPartialWrtSunGravitationalParameter = calculateAccelerationWrtParameterPartials(
                gravitationalParameterParameter, gravitationalAcceleration, 1.0E16 );
    Eigen::Vector3d testPartialWrtEarthGravitationalParameter = calculateAccelerationWrtParameterPartials(
                earthGravitationalParameterParameter, gravitationalAcceleration, 1.0E16 );
    Eigen::Vector3d testPartialWrtMoonGravitationalParameter = calculateAccelerationWrtParameterPartials(
                moonGravitationalParameterParameter, gravitationalAcceleration, 1.0E16 );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonPosition,
                                       partialWrtMoonPosition, 1.0E-7 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonVelocity,
                                       partialWrtMoonVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                       partialWrtSunPosition, 1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                       partialWrtSunVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                       partialWrtEarthPosition, 1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                       partialWrtEarthVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunGravitationalParameter,
                                       partialWrtSunGravitationalParameter, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonGravitationalParameter,
                                       partialWrtMoonGravitationalParameter, std::numeric_limits< double >::epsilon(  ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthGravitationalParameter,
                                       partialWrtEarthGravitationalParameter, std::numeric_limits< double >::epsilon(  ) );
}


void updateFlightConditionsWithPerturbedState(
        const std::shared_ptr< aerodynamics::FlightConditions > flightConditions,
        const double timeToUpdate )
{
    flightConditions->resetCurrentTime( TUDAT_NAN );
    flightConditions->updateConditions( timeToUpdate );
}

BOOST_AUTO_TEST_CASE( testAerodynamicAccelerationPartials )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    using namespace tudat;
    // Create Earth object
    std::map< std::string, std::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings( { "Earth" } );
    defaultBodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ) );
    NamedBodyMap bodyMap = createBodies( defaultBodySettings );

    // Create vehicle objects.
    double vehicleMass = 5.0E3;
    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );


    bool areCoefficientsInAerodynamicFrame = 1;
    Eigen::Vector3d aerodynamicCoefficients = ( Eigen::Vector3d( ) << 2.5, -0.1, 0.5 ).finished( );

    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                2.0, 4.0, 1.5, Eigen::Vector3d::Zero( ), aerodynamicCoefficients, Eigen::Vector3d::Zero( ),
                areCoefficientsInAerodynamicFrame, 1 );
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Set spherical elements for vehicle.
    Eigen::Vector6d vehicleSphericalEntryState;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.0;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.7E3;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            -0.9 * mathematical_constants::PI / 180.0;
    vehicleSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

    // Convert vehicle state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                vehicleSphericalEntryState );

    bodyMap.at( "Earth" )->setStateFromEphemeris( 0.0 );
    bodyMap.at( "Vehicle" )->setState( systemInitialState );


    std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel =
            simulation_setup::createAerodynamicAcceleratioModel(
                bodyMap[ "Vehicle" ], bodyMap[ "Earth" ], "Vehicle", "Earth" );
    bodyMap.at( "Vehicle" )->getFlightConditions( )->updateConditions( 0.0 );
    accelerationModel->updateMembers( 0.0 );

    std::shared_ptr< AccelerationPartial > aerodynamicAccelerationPartial =
            createAnalyticalAccelerationPartial(
                accelerationModel, std::make_pair( "Vehicle", bodyMap[ "Vehicle" ] ),
            std::make_pair( "Earth", bodyMap[ "Earth" ] ), bodyMap );

    // Create gravitational parameter object.
    std::shared_ptr< EstimatableParameter< double > > dragCoefficientParameter = std::make_shared<
            ConstantDragCoefficient >( std::dynamic_pointer_cast< aerodynamics::CustomAerodynamicCoefficientInterface >(
                                           bodyMap[ "Vehicle" ]->getAerodynamicCoefficientInterface( ) ), "Vehicle" );

    // Calculate analytical partials.
    aerodynamicAccelerationPartial->update( 0.0 );
    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    aerodynamicAccelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    aerodynamicAccelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    aerodynamicAccelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    aerodynamicAccelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::Vector3d partialWrtDragCoefficient = aerodynamicAccelerationPartial->wrtParameter(
                dragCoefficientParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );

    std::function< void( ) > environmentUpdateFunction =
            std::bind( &updateFlightConditionsWithPerturbedState, bodyMap.at( "Vehicle" )->getFlightConditions( ), 0.0 );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1.0, 1.0, 1.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            std::bind( &Body::setState, bodyMap.at( "Vehicle" ), std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, bodyMap.at( "Earth" ), std::placeholders::_1 );
    std::function< Eigen::Vector6d ( ) > vehicleStateGetFunction =
            std::bind( &Body::getState, bodyMap.at( "Vehicle" ) );
    std::function< Eigen::Vector6d ( ) > earthStateGetFunction =
            std::bind( &Body::getState, bodyMap.at( "Earth" ) );

    // Calculate numerical partials.
    testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, bodyMap.at( "Vehicle" )->getState( ), positionPerturbation, 0,
                environmentUpdateFunction);
    testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, bodyMap.at( "Vehicle" )->getState( ), velocityPerturbation, 3,
                environmentUpdateFunction );
    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, accelerationModel, bodyMap.at( "Earth" )->getState( ), positionPerturbation, 0,
                environmentUpdateFunction );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, accelerationModel, bodyMap.at( "Earth" )->getState( ), velocityPerturbation, 3,
                environmentUpdateFunction );

    Eigen::Vector3d testPartialWrtDragCoefficient = calculateAccelerationWrtParameterPartials(
                dragCoefficientParameter, accelerationModel, 1.0E-4, environmentUpdateFunction );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                       partialWrtVehiclePosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                       partialWrtVehicleVelocity, 1.0E-6  );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                       partialWrtEarthPosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                       partialWrtEarthVelocity, 1.0E-6 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtDragCoefficient,
                                       partialWrtDragCoefficient, 1.0E-10 );
}



BOOST_AUTO_TEST_CASE( testIndividualTermsOfRelativisticAccelerationPartial )
{
    // Create earth and vehicle bodies.
    std::shared_ptr< Body > sun = std::make_shared< Body >( );
    std::shared_ptr< Body > mercury = std::make_shared< Body >( );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            std::bind( &Body::setState, sun, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > mercuryStateSetFunction =
            std::bind( &Body::setState, mercury, std::placeholders::_1  );
    std::function< Eigen::Vector6d( ) > sunStateGetFunction =
            std::bind( &Body::getState, sun );
    std::function< Eigen::Vector6d( ) > mercuryStateGetFunction =
            std::bind( &Body::getState, mercury);

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Set vehicle and earth state.
    sun->setState( getBodyCartesianStateAtEpoch(  "Sun", "SSB", "J2000", "NONE", 1.0E6 ) );
    mercury->setState( getBodyCartesianStateAtEpoch(  "Mercury", "SSB", "J2000", "NONE", 1.0E6 ) );

    NamedBodyMap bodyMap;
    bodyMap[ "Sun" ] = sun;
    bodyMap[ "Mercury" ] = mercury;

    // Create gravity field.
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );
    std::shared_ptr< gravitation::GravityFieldModel > sunGravityField =
            createGravityFieldModel( gravityFieldSettings, "Sun", bodyMap );
    sun->setGravityFieldModel( sunGravityField );
    std::shared_ptr< gravitation::GravityFieldModel > mercuryGravityField =
            createGravityFieldModel( gravityFieldSettings, "Mercury", bodyMap );
    sun->setGravityFieldModel( mercuryGravityField );

    // Get angular momentum Sun
    const Eigen::Vector3d sunAngularMomentumVector(0.0, 0.0, 190.0E39);
    const Eigen::Vector3d sunAngularMomentumVectorPerUnitMass =
            sunAngularMomentumVector /
            ( sunGravityField->getGravitationalParameter() / physical_constants::GRAVITATIONAL_CONSTANT );

    std::function< Eigen::Vector3d( ) > angularMomentumFunction;
    angularMomentumFunction = [ = ]( ){ return sunAngularMomentumVectorPerUnitMass; };

    std::function< std::string( ) > nameOfBodyExertingAccelerationFunction;
    nameOfBodyExertingAccelerationFunction = [ = ]( ){ return "Sun"; };

    // Create acceleration model.
    std::function< double( ) > ppnParameterGammaFunction = std::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet );
    std::function< double( ) > ppnParameterBetaFunction = std::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet );
    std::function< double( ) > ppnParameterAlpha1Function = std::bind( &PPNParameterSet::getParameterAlpha1, ppnParameterSet );
    std::function< double( ) > ppnParameterAlpha2Function = std::bind( &PPNParameterSet::getParameterAlpha2, ppnParameterSet );

    std::shared_ptr< RelativisticAccelerationCorrection > accelerationModel =
            std::make_shared< RelativisticAccelerationCorrection >
            ( std::bind( &Body::getState, mercury ),
              std::bind( &Body::getState, sun ),
              std::bind( &GravityFieldModel::getGravitationalParameter, sunGravityField ),
              std::bind( &GravityFieldModel::getGravitationalParameter, mercuryGravityField ),
              angularMomentumFunction, nameOfBodyExertingAccelerationFunction,
              ppnParameterGammaFunction, ppnParameterBetaFunction,
              ppnParameterAlpha1Function, ppnParameterAlpha2Function);

    std::cout<<"ss:"<<accelerationModel->getCalculateSchwarzschildCorrection()
             <<" lt:"<<accelerationModel->getCalculateLenseThirringCorrection()
             <<" ds:"<<accelerationModel->getCalculateDeSitterCorrection()<<std::endl;

    // Create acceleration partial object.
    std::shared_ptr< RelativisticAccelerationPartial > accelerationPartial = std::make_shared< RelativisticAccelerationPartial >(
                accelerationModel, "Mercury", "Sun" );

    // Create parameter objects.
    std::shared_ptr< EstimatableParameter< double > > gravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( sunGravityField, "Sun" );
    std::shared_ptr< EstimatableParameter< double > > ppnParameterGamma = std::make_shared<
            PPNParameterGamma >( ppnParameterSet );
    std::shared_ptr< EstimatableParameter< double > > ppnParameterBeta = std::make_shared<
            PPNParameterBeta >( ppnParameterSet );
    std::shared_ptr< EstimatableParameter< double > > ppnParameterAlpha1 = std::make_shared<
            PPNParameterAlpha1 >( ppnParameterSet );
    std::shared_ptr< EstimatableParameter< double > > ppnParameterAlpha2 = std::make_shared<
            PPNParameterAlpha2 >( ppnParameterSet );

    // Calculate analytical partials.
    accelerationModel->updateMembers( );
    accelerationPartial->update( );

    std::vector< Eigen::MatrixXd > partialsWrtSunPosition;
    std::vector< Eigen::MatrixXd > partialsWrtSunVelocity;
    std::vector< Eigen::MatrixXd > partialsWrtMercuryPosition;
    std::vector< Eigen::MatrixXd > partialsWrtMercuryVelocity;

    for (unsigned int i=1; i<4; i++){
        partialsWrtSunPosition.push_back(
                    -1.0*accelerationPartial->returnPositionPartialOfOneTerm(i));
        partialsWrtSunVelocity.push_back(
                    -1.0*accelerationPartial->returnVelocityPartialOfOneTerm(i));
        partialsWrtMercuryPosition.push_back(
                    accelerationPartial->returnPositionPartialOfOneTerm(i));
        partialsWrtMercuryVelocity.push_back(
                    accelerationPartial->returnVelocityPartialOfOneTerm(i));
    }

    Eigen::Vector3d partialWrtSunGravitationalParameter = accelerationPartial->wrtParameter(
                gravitationalParameterParameter );

    Eigen::Vector3d partialWrtGamma = accelerationPartial->wrtParameter( ppnParameterGamma );
    Eigen::Vector3d partialWrtBeta = accelerationPartial->wrtParameter( ppnParameterBeta );
    Eigen::Vector3d partialWrtAlpha1 = accelerationPartial->wrtParameter( ppnParameterAlpha1 );
    Eigen::Vector3d partialWrtAlpha2 = accelerationPartial->wrtParameter( ppnParameterAlpha2 );

    // Declare perturbations in position for numerical partial
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10.0, 10.0, 10.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;


    // Calculate numerical partials.
    std::vector< Eigen::MatrixXd > testPartialsWrtSunPosition;
    std::vector< Eigen::MatrixXd > testPartialsWrtSunVelocity;
    std::vector< Eigen::MatrixXd > testPartialsWrtMercuryPosition;
    std::vector< Eigen::MatrixXd > testPartialsWrtMercuryVelocity;

    for (unsigned int i=1; i<4; i++){
        testPartialsWrtSunPosition.push_back(
                    calculateAccelerationWrtStatePartialsRelativity(
                                    sunStateSetFunction, accelerationModel, sun->getState( ), positionPerturbation, 0, i ));
        testPartialsWrtSunVelocity.push_back(
                    calculateAccelerationWrtStatePartialsRelativity(
                                    sunStateSetFunction, accelerationModel, sun->getState( ), velocityPerturbation, 3, i ));
        testPartialsWrtMercuryPosition.push_back(
                    calculateAccelerationWrtStatePartialsRelativity(
                                    mercuryStateSetFunction, accelerationModel, mercury->getState( ), positionPerturbation, 0, i ));
        testPartialsWrtMercuryVelocity.push_back(
                    calculateAccelerationWrtStatePartialsRelativity(
                                    mercuryStateSetFunction, accelerationModel, mercury->getState( ), positionPerturbation, 3, i ));
    }

    Eigen::Vector3d testPartialWrtSunGravitationalParameter = calculateAccelerationWrtParameterPartials(
                gravitationalParameterParameter, accelerationModel, 1.0E10 );

    Eigen::Vector3d testPartialWrtPpnParameterGamma = calculateAccelerationWrtParameterPartials(
                ppnParameterGamma, accelerationModel, 100.0 );
    Eigen::Vector3d testPartialWrtPpnParameterBeta = calculateAccelerationWrtParameterPartials(
                ppnParameterBeta, accelerationModel, 100.0 );

    Eigen::Vector3d testPartialWrtPpnParameterAlpha1 = calculateAccelerationWrtParameterPartials(
                ppnParameterAlpha1, accelerationModel, 100.0 );
    Eigen::Vector3d testPartialWrtPpnParameterAlpha2 = calculateAccelerationWrtParameterPartials(
                ppnParameterAlpha2, accelerationModel, 100.0 );


    double stateTolerance = 1.0E-3; //orginally E-7?
    double parameterTolerance = 1.0E-4; //originally E-8

    // Compare numerical and analytical results.
    for (unsigned int i=1; i<4; i++){

        switch(i){
            case 1: {std::cout<<"---- Conventional Schwarzschild correction:"<<std::endl; break;}
            case 2: {std::cout<<"---- Schwarzschild alpha correction:"<<std::endl; break;}
            case 3: {std::cout<<"---- Lense-thirring correction:"<<std::endl; break;}
        }

        std::cout<<"wrt sun position:"
                 <<std::endl<<partialsWrtSunPosition.at(i-1)
                 <<std::endl<<testPartialsWrtSunPosition.at(i-1)
                 <<std::endl<<partialsWrtSunPosition.at(i-1)-testPartialsWrtSunPosition.at(i-1)
                 <<std::endl;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialsWrtSunPosition.at(i-1),
                                           partialsWrtSunPosition.at(i-1), stateTolerance );

        std::cout<<"wrt sun velocity:"
                 <<std::endl<<partialsWrtSunVelocity.at(i-1)
                 <<std::endl<<testPartialsWrtSunVelocity.at(i-1)
                 <<std::endl<<partialsWrtSunVelocity.at(i-1)-testPartialsWrtSunVelocity.at(i-1)
                 <<std::endl;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialsWrtSunVelocity.at(i-1),
                                           partialsWrtSunVelocity.at(i-1), stateTolerance );

        std::cout<<"wrt mercury position:"
                 <<std::endl<<partialsWrtMercuryPosition.at(i-1)
                 <<std::endl<<testPartialsWrtMercuryPosition.at(i-1)
                 <<std::endl<<partialsWrtMercuryPosition.at(i-1)-testPartialsWrtMercuryPosition.at(i-1)
                 <<std::endl;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialsWrtMercuryPosition.at(i-1),
                                           partialsWrtMercuryPosition.at(i-1), stateTolerance );

        std::cout<<"wrt mercury velocity:"
                 <<std::endl<<partialsWrtMercuryVelocity.at(i-1)
                 <<std::endl<<testPartialsWrtMercuryVelocity.at(i-1)
                 <<std::endl<<partialsWrtMercuryVelocity.at(i-1)-testPartialsWrtMercuryVelocity.at(i-1)
                 <<std::endl;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialsWrtMercuryVelocity.at(i-1),
                                           partialsWrtMercuryVelocity.at(i-1), stateTolerance );


    }

    std::cout<<"wrt grav parameter:"
             <<std::endl<<partialWrtSunGravitationalParameter.transpose()
             <<std::endl<<testPartialWrtSunGravitationalParameter.transpose()
             <<std::endl<<partialWrtSunGravitationalParameter.transpose()-testPartialWrtSunGravitationalParameter.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunGravitationalParameter,
                                       partialWrtSunGravitationalParameter, parameterTolerance );

    std::cout<<"wrt gamma:"
             <<std::endl<<partialWrtGamma.transpose()
             <<std::endl<<testPartialWrtPpnParameterGamma.transpose()
             <<std::endl<<partialWrtGamma.transpose()-testPartialWrtPpnParameterGamma.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterGamma, partialWrtGamma, parameterTolerance );

    std::cout<<"wrt beta:"
             <<std::endl<<partialWrtBeta.transpose()
             <<std::endl<<testPartialWrtPpnParameterBeta.transpose()
             <<std::endl<<partialWrtBeta.transpose()-testPartialWrtPpnParameterBeta.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterBeta, partialWrtBeta, parameterTolerance );

    std::cout<<"wrt alpha1:"
             <<std::endl<<partialWrtAlpha1.transpose()
             <<std::endl<<testPartialWrtPpnParameterAlpha1.transpose()
             <<std::endl<<partialWrtAlpha1.transpose()-testPartialWrtPpnParameterAlpha1.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterAlpha1, partialWrtAlpha1, parameterTolerance );

    std::cout<<"wrt alpha2:"
             <<std::endl<<partialWrtAlpha2.transpose()
             <<std::endl<<testPartialWrtPpnParameterAlpha2.transpose()
             <<std::endl<<partialWrtAlpha2.transpose()-testPartialWrtPpnParameterAlpha2.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterAlpha2, partialWrtAlpha2, parameterTolerance );
}





BOOST_AUTO_TEST_CASE( testRelativisticAccelerationPartial )
{
    // Create earth and vehicle bodies.
    std::shared_ptr< Body > sun = std::make_shared< Body >( );
    std::shared_ptr< Body > mercury = std::make_shared< Body >( );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            std::bind( &Body::setState, sun, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > mercuryStateSetFunction =
            std::bind( &Body::setState, mercury, std::placeholders::_1  );
    std::function< Eigen::Vector6d( ) > sunStateGetFunction =
            std::bind( &Body::getState, sun );
    std::function< Eigen::Vector6d( ) > mercuryStateGetFunction =
            std::bind( &Body::getState, mercury);

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Set vehicle and earth state.
    sun->setState( getBodyCartesianStateAtEpoch(  "Sun", "SSB", "J2000", "NONE", 1.0E6 ) );
    mercury->setState( getBodyCartesianStateAtEpoch(  "Mercury", "SSB", "J2000", "NONE", 1.0E6 ) );

    NamedBodyMap bodyMap;
    bodyMap[ "Sun" ] = sun;
    bodyMap[ "Mercury" ] = mercury;

    // Create gravity field.
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );
    std::shared_ptr< gravitation::GravityFieldModel > sunGravityField =
            createGravityFieldModel( gravityFieldSettings, "Sun", bodyMap );
    sun->setGravityFieldModel( sunGravityField );
    std::shared_ptr< gravitation::GravityFieldModel > mercuryGravityField =
            createGravityFieldModel( gravityFieldSettings, "Mercury", bodyMap );
    sun->setGravityFieldModel( mercuryGravityField );

    // Get angular momentum Sun
    const Eigen::Vector3d sunAngularMomentumVector(0.0, 0.0, 190.0E39);
    const Eigen::Vector3d sunAngularMomentumVectorPerUnitMass =
            sunAngularMomentumVector /
            ( sunGravityField->getGravitationalParameter() / physical_constants::GRAVITATIONAL_CONSTANT );

    std::function< Eigen::Vector3d( ) > angularMomentumFunction;
    angularMomentumFunction = [ = ]( ){ return sunAngularMomentumVectorPerUnitMass; };

    std::function< std::string( ) > nameOfBodyExertingAccelerationFunction;
    nameOfBodyExertingAccelerationFunction = [ = ]( ){ return "Sun"; };

    // Create acceleration model.
    std::function< double( ) > ppnParameterGammaFunction = std::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet );
    std::function< double( ) > ppnParameterBetaFunction = std::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet );
    std::function< double( ) > ppnParameterAlpha1Function = std::bind( &PPNParameterSet::getParameterAlpha1, ppnParameterSet );
    std::function< double( ) > ppnParameterAlpha2Function = std::bind( &PPNParameterSet::getParameterAlpha2, ppnParameterSet );

    std::shared_ptr< RelativisticAccelerationCorrection > accelerationModel =
            std::make_shared< RelativisticAccelerationCorrection >
            ( std::bind( &Body::getState, mercury ),
              std::bind( &Body::getState, sun ),
              std::bind( &GravityFieldModel::getGravitationalParameter, sunGravityField ),
              std::bind( &GravityFieldModel::getGravitationalParameter, mercuryGravityField ),
              angularMomentumFunction, nameOfBodyExertingAccelerationFunction,
              ppnParameterGammaFunction, ppnParameterBetaFunction,
              ppnParameterAlpha1Function, ppnParameterAlpha2Function);

    std::cout<<"ss:"<<accelerationModel->getCalculateSchwarzschildCorrection()
             <<" lt:"<<accelerationModel->getCalculateLenseThirringCorrection()
             <<" ds:"<<accelerationModel->getCalculateDeSitterCorrection()<<std::endl;

    // Create acceleration partial object.
    std::shared_ptr< RelativisticAccelerationPartial > accelerationPartial = std::make_shared< RelativisticAccelerationPartial >(
                accelerationModel, "Mercury", "Sun" );

    // Create parameter objects.
    std::shared_ptr< EstimatableParameter< double > > gravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( sunGravityField, "Sun" );
    std::shared_ptr< EstimatableParameter< double > > ppnParameterGamma = std::make_shared<
            PPNParameterGamma >( ppnParameterSet );
    std::shared_ptr< EstimatableParameter< double > > ppnParameterBeta = std::make_shared<
            PPNParameterBeta >( ppnParameterSet );
    std::shared_ptr< EstimatableParameter< double > > ppnParameterAlpha1 = std::make_shared<
            PPNParameterAlpha1 >( ppnParameterSet );
    std::shared_ptr< EstimatableParameter< double > > ppnParameterAlpha2 = std::make_shared<
            PPNParameterAlpha2 >( ppnParameterSet );

    // Calculate analytical partials.
    accelerationModel->updateMembers( );
    accelerationPartial->update( );

    Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtMercuryPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtMercuryPosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtMercuryVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtMercuryVelocity.block( 0, 0, 3, 3 ) );

    Eigen::Vector3d partialWrtSunGravitationalParameter = accelerationPartial->wrtParameter(
                gravitationalParameterParameter );

    Eigen::Vector3d partialWrtGamma = accelerationPartial->wrtParameter( ppnParameterGamma );
    Eigen::Vector3d partialWrtBeta = accelerationPartial->wrtParameter( ppnParameterBeta );
    Eigen::Vector3d partialWrtAlpha1 = accelerationPartial->wrtParameter( ppnParameterAlpha1 );
    Eigen::Vector3d partialWrtAlpha2 = accelerationPartial->wrtParameter( ppnParameterAlpha2 );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10.0, 10.0, 10.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Calculate numerical partials.
    Eigen::Matrix3d testPartialWrtMercuryPosition = calculateAccelerationWrtStatePartials(
                mercuryStateSetFunction, accelerationModel, mercury->getState( ), positionPerturbation, 0 );
    Eigen::Matrix3d testPartialWrtMercuryVelocity = calculateAccelerationWrtStatePartials(
                mercuryStateSetFunction, accelerationModel, mercury->getState( ), velocityPerturbation, 3 );

    Eigen::Matrix3d testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, accelerationModel, sun->getState( ), positionPerturbation, 0 );
    Eigen::Matrix3d testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, accelerationModel, sun->getState( ), velocityPerturbation, 3 );

    Eigen::Vector3d testPartialWrtSunGravitationalParameter = calculateAccelerationWrtParameterPartials(
                gravitationalParameterParameter, accelerationModel, 1.0E10 );

    Eigen::Vector3d testPartialWrtPpnParameterGamma = calculateAccelerationWrtParameterPartials(
                ppnParameterGamma, accelerationModel, 100.0 );
    Eigen::Vector3d testPartialWrtPpnParameterBeta = calculateAccelerationWrtParameterPartials(
                ppnParameterBeta, accelerationModel, 100.0 );

    Eigen::Vector3d testPartialWrtPpnParameterAlpha1 = calculateAccelerationWrtParameterPartials(
                ppnParameterAlpha1, accelerationModel, 100.0 );
    Eigen::Vector3d testPartialWrtPpnParameterAlpha2 = calculateAccelerationWrtParameterPartials(
                ppnParameterAlpha2, accelerationModel, 100.0 );



    // Compare numerical and analytical results.
    double stateTolerance = 1.0E-3; //orginally E-7?
    double parameterTolerance = 1.0E-4; //originally E-8


    std::cout<<"wrt sun position:"
             <<std::endl<<partialWrtSunPosition
             <<std::endl<<testPartialWrtSunPosition
             <<std::endl<<partialWrtSunPosition-testPartialWrtSunPosition
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                       partialWrtSunPosition, stateTolerance);

    std::cout<<"wrt sun velocity:"
             <<std::endl<<partialWrtSunVelocity
             <<std::endl<<testPartialWrtSunVelocity
             <<std::endl<<partialWrtSunVelocity-testPartialWrtSunVelocity
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                       partialWrtSunVelocity, stateTolerance);

    std::cout<<"wrt mercury position:"
             <<std::endl<<partialWrtMercuryPosition
             <<std::endl<<testPartialWrtMercuryPosition
             <<std::endl<<partialWrtMercuryPosition-testPartialWrtMercuryPosition
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMercuryPosition,
                                       partialWrtMercuryPosition, stateTolerance );

    std::cout<<"wrt mercury velocity:"
             <<std::endl<<partialWrtMercuryVelocity
             <<std::endl<<testPartialWrtMercuryVelocity
             <<std::endl<<partialWrtMercuryVelocity-testPartialWrtMercuryVelocity
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMercuryVelocity,
                                       partialWrtMercuryVelocity, stateTolerance );

    std::cout<<"wrt grav parameter:"
             <<std::endl<<partialWrtSunGravitationalParameter.transpose()
             <<std::endl<<testPartialWrtSunGravitationalParameter.transpose()
             <<std::endl<<partialWrtSunGravitationalParameter.transpose()-testPartialWrtSunGravitationalParameter.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunGravitationalParameter,
                                       partialWrtSunGravitationalParameter, parameterTolerance );

    std::cout<<"wrt gamma:"
             <<std::endl<<partialWrtGamma.transpose()
             <<std::endl<<testPartialWrtPpnParameterGamma.transpose()
             <<std::endl<<partialWrtGamma.transpose()-testPartialWrtPpnParameterGamma.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterGamma, partialWrtGamma, parameterTolerance );

    std::cout<<"wrt beta:"
             <<std::endl<<partialWrtBeta.transpose()
             <<std::endl<<testPartialWrtPpnParameterBeta.transpose()
             <<std::endl<<partialWrtBeta.transpose()-testPartialWrtPpnParameterBeta.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterBeta, partialWrtBeta, parameterTolerance );

    std::cout<<"wrt alpha1:"
             <<std::endl<<partialWrtAlpha1.transpose()
             <<std::endl<<testPartialWrtPpnParameterAlpha1.transpose()
             <<std::endl<<partialWrtAlpha1.transpose()-testPartialWrtPpnParameterAlpha1.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterAlpha1, partialWrtAlpha1, parameterTolerance );

    std::cout<<"wrt alpha2:"
             <<std::endl<<partialWrtAlpha2.transpose()
             <<std::endl<<testPartialWrtPpnParameterAlpha2.transpose()
             <<std::endl<<partialWrtAlpha2.transpose()-testPartialWrtPpnParameterAlpha2.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterAlpha2, partialWrtAlpha2, parameterTolerance );
}






BOOST_AUTO_TEST_CASE( testSEPViolationAccelerationPartial )
{

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );

    // Define bodies
    unsigned int totalNumberOfBodies = 8;
    std::vector< std::string > bodyNames;
    bodyNames.resize( totalNumberOfBodies );
    bodyNames[ 0 ] = "Sun";
    bodyNames[ 1 ] = "Mercury";
    bodyNames[ 2 ] = "Venus";
    bodyNames[ 3 ] = "Earth";
    bodyNames[ 4 ] = "Mars";
    bodyNames[ 5 ] = "Jupiter";
    bodyNames[ 6 ] = "Saturn";
//    bodyNames[ 7 ] = "Uranus";
//    bodyNames[ 8 ] = "Neptune";
    bodyNames[ 7 ] = "Moon";

    // Manually load body settings
    NamedBodyMap bodyMap;
    for (unsigned int i=0; i<totalNumberOfBodies; i++){

        std::shared_ptr< Body > body = std::make_shared< Body >( );

        std::function< void( Eigen::Vector6d ) > bodyStateSetFunction =
                std::bind( &Body::setState, body, std::placeholders::_1  );
        std::function< Eigen::Vector6d( ) > bodyStateGetFunction =
                std::bind( &Body::getState, body );

        body->setState( getBodyCartesianStateAtEpoch( bodyNames.at(i), "SSB", "J2000", "NONE", 1.0E6 ) );

        bodyMap[ bodyNames.at(i) ] = body;

        std::shared_ptr< gravitation::GravityFieldModel > bodyGravityField =
                createGravityFieldModel( gravityFieldSettings, bodyNames.at(i), bodyMap );
        body->setGravityFieldModel( bodyGravityField );

        bodyMap[ bodyNames.at(i) ] = body;

    }




//    // Create earth and vehicle bodies.
//    std::shared_ptr< Body > sun = std::make_shared< Body >( );
//    std::shared_ptr< Body > mercury = std::make_shared< Body >( );
//    std::shared_ptr< Body > venus = std::make_shared< Body >( );
//    std::shared_ptr< Body > jupiter = std::make_shared< Body >( );

//    // Create links to set and get state functions of bodies.
//    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
//            std::bind( &Body::setState, sun, std::placeholders::_1  );
//    std::function< void( Eigen::Vector6d ) > mercuryStateSetFunction =
//            std::bind( &Body::setState, mercury, std::placeholders::_1  );
//    std::function< void( Eigen::Vector6d ) > venusStateSetFunction =
//            std::bind( &Body::setState, venus, std::placeholders::_1  );
//    std::function< void( Eigen::Vector6d ) > jupiterStateSetFunction =
//            std::bind( &Body::setState, jupiter, std::placeholders::_1  );

//    std::function< Eigen::Vector6d( ) > sunStateGetFunction =
//            std::bind( &Body::getState, sun );
//    std::function< Eigen::Vector6d( ) > mercuryStateGetFunction =
//            std::bind( &Body::getState, mercury);
//    std::function< Eigen::Vector6d( ) > venusStateGetFunction =
//            std::bind( &Body::getState, venus );
//    std::function< Eigen::Vector6d( ) > jupiterStateGetFunction =
//            std::bind( &Body::getState, jupiter);

//    // Set vehicle and earth state.
//    sun->setState( getBodyCartesianStateAtEpoch(  "Sun", "SSB", "J2000", "NONE", 1.0E6 ) );
//    mercury->setState( getBodyCartesianStateAtEpoch(  "Mercury", "SSB", "J2000", "NONE", 1.0E6 ) );
//    venus->setState( getBodyCartesianStateAtEpoch(  "Venus", "SSB", "J2000", "NONE", 1.0E6 ) );
//    jupiter->setState( getBodyCartesianStateAtEpoch(  "Jupiter", "SSB", "J2000", "NONE", 1.0E6 ) );

//    NamedBodyMap bodyMap;
//    bodyMap[ "Sun" ] = sun;
//    bodyMap[ "Mercury" ] = mercury;
//    bodyMap[ "Venus" ] = venus;
//    bodyMap[ "Jupiter" ] = jupiter;

//    std::vector< std::string > bodyNames;
//    bodyNames.push_back( "Sun" );
//    bodyNames.push_back( "Mercury" );
//    bodyNames.push_back( "Venus" );
//    bodyNames.push_back( "Jupiter" );

//    // Create gravity field.
//    std::shared_ptr< gravitation::GravityFieldModel > sunGravityField =
//            createGravityFieldModel( gravityFieldSettings, "Sun", bodyMap );
//    sun->setGravityFieldModel( sunGravityField );
//    std::shared_ptr< gravitation::GravityFieldModel > mercuryGravityField =
//            createGravityFieldModel( gravityFieldSettings, "Mercury", bodyMap );
//    mercury->setGravityFieldModel( mercuryGravityField );
//    std::shared_ptr< gravitation::GravityFieldModel > venusGravityField =
//            createGravityFieldModel( gravityFieldSettings, "Venus", bodyMap );
//    venus->setGravityFieldModel( venusGravityField );
//    std::shared_ptr< gravitation::GravityFieldModel > jupiterGravityField =
//            createGravityFieldModel( gravityFieldSettings, "Jupiter", bodyMap );
//    jupiter->setGravityFieldModel( jupiterGravityField );

    std::function< Eigen::Vector3d( ) > nordtvedtPartialFunction =
            std::bind( getNordtvedtPartial,
                       bodyMap.at( "Sun" ), bodyMap.at( "Mercury" ),
                       "Sun", "Mercury",
                       bodyMap, bodyNames);

    // run unit test for two cases: nordtvedt constraint true and false
    for (int j = 0; j<2; j++){

        std::function< bool( ) > useNordtvedtConstraintFunction;
        std::function< double( ) >  nordtvedtParameterFunction;
        if(j){
            useNordtvedtConstraintFunction = [ = ]( ){ return true; };
//            ppnParameterSet->setParameterGamma(1.0);
//            ppnParameterSet->setParameterBeta(1.0);
//            ppnParameterSet->setParameterAlpha1(0.0);
//            ppnParameterSet->setParameterAlpha2(0.0);
            ppnParameterSet->setParameterGamma(1.0 + 1.0E-4);
            ppnParameterSet->setParameterBeta(1.0 + 2.0E-4);
            ppnParameterSet->setParameterAlpha1(3.0E-4);
            ppnParameterSet->setParameterAlpha2(4.0E-4);
            nordtvedtParameterFunction =
                    std::bind( &PPNParameterSet::getNordtvedtParameterFromPpnParameters, ppnParameterSet );

        } else{
            useNordtvedtConstraintFunction = [ = ]( ){ return false; };
//            ppnParameterSet->setNordtvedtParameter(0.0);
            ppnParameterSet->setNordtvedtParameter(1.0E-3);
            nordtvedtParameterFunction =
                    std::bind( &PPNParameterSet::getNordtvedtParameter, ppnParameterSet );

        }

        std::cout<<"nordtvedt constraint true/false: "<<useNordtvedtConstraintFunction()<<std::endl;

        std::function< Eigen::Vector3d( ) > sepPositionCorrectionFunction =
                std::bind( getSEPCorrectedPosition,
                           bodyMap.at( "Sun" ), "Sun",
                           bodyMap, bodyNames,
                           nordtvedtParameterFunction);


        std::function< Eigen::Vector3d( ) > sepPositionCorrectionFunctionUnsimplified =
                std::bind( getSEPCorrectedPositionUnsimplified,
                           bodyMap.at( "Sun" ), "Sun",
                           bodyMap, bodyNames,
                           nordtvedtParameterFunction);

//        std::cout<<useNordtvedtConstraintFunction()<<std::endl;
//        std::cout<<sepPositionCorrectionFunction().transpose()<<std::endl;
//        std::cout<<nordtvedtPartialFunction().transpose()<<std::endl;

        // Create acceleration model.
        std::shared_ptr< SEPViolationAcceleration > accelerationModel
                = std::make_shared< SEPViolationAcceleration >
                ( std::bind( &Body::getPosition, bodyMap.at( "Mercury" ) ),
                  std::bind( &Body::getPosition, bodyMap.at( "Sun" ) ),
                  sepPositionCorrectionFunction,
                  std::bind( &GravityFieldModel::getGravitationalParameter, bodyMap.at( "Sun" )->getGravityFieldModel() ),
                  nordtvedtPartialFunction,
                  useNordtvedtConstraintFunction,
                  nordtvedtParameterFunction
                  );

        // Create acceleration model.
        std::shared_ptr< SEPViolationAcceleration > accelerationModelUnsimplified
                = std::make_shared< SEPViolationAcceleration >
                ( std::bind( &Body::getPosition, bodyMap.at( "Mercury" ) ),
                  std::bind( &Body::getPosition, bodyMap.at( "Sun" ) ),
                  sepPositionCorrectionFunctionUnsimplified,
                  std::bind( &GravityFieldModel::getGravitationalParameter, bodyMap.at( "Sun" )->getGravityFieldModel() ),
                  nordtvedtPartialFunction,
                  useNordtvedtConstraintFunction,
                  nordtvedtParameterFunction
                  );


        // Create acceleration partial object.
        std::shared_ptr< SEPViolationAccelerationPartial > accelerationPartial =
                std::make_shared< SEPViolationAccelerationPartial >(
                    accelerationModel, "Mercury", "Sun" );


        // Create parameter objects.
        std::shared_ptr< EstimatableParameter< double > > gravitationalParameterParameter = std::make_shared<
                GravitationalParameter >( bodyMap.at( "Sun" )->getGravityFieldModel(), "Sun" );
        std::shared_ptr< EstimatableParameter< double > > ppnParameterGamma = std::make_shared<
                PPNParameterGamma >( ppnParameterSet );
        std::shared_ptr< EstimatableParameter< double > > ppnParameterBeta = std::make_shared<
                PPNParameterBeta >( ppnParameterSet );
        std::shared_ptr< EstimatableParameter< double > > ppnParameterAlpha1 = std::make_shared<
                PPNParameterAlpha1 >( ppnParameterSet );
        std::shared_ptr< EstimatableParameter< double > > ppnParameterAlpha2 = std::make_shared<
                PPNParameterAlpha2 >( ppnParameterSet );
        std::shared_ptr< EstimatableParameter< double > > nordtvedtParameter = std::make_shared<
                PPNNordtvedtParameter >( ppnParameterSet );


        // Calculate analytical partials.
        accelerationModel->updateMembers( );
        accelerationModelUnsimplified->updateMembers( );
        accelerationPartial->update( );

//        Eigen::Vector3d currentSEPcorrection = accelerationModel->;
//        Eigen::Vector3d currentSEPcorrectionUnsimplified = accelerationModelUnsimplified->getSEPPositionCorrectionFunction();

        Eigen::Vector3d currentAcceleration = accelerationModel->getAcceleration();
        Eigen::Vector3d currentAccelerationUnsimplified = accelerationModelUnsimplified->getAcceleration();

        std::cout<<"acceleration: "<<currentAcceleration.transpose()<<std::endl;
        std::cout<<"acceleration unsimplified: "<<currentAccelerationUnsimplified.transpose()<<std::endl;
        std::cout<<"difference: "<<(currentAcceleration-currentAccelerationUnsimplified).transpose()<<std::endl;
        std::cout<<"relative difference: "<<((currentAcceleration-currentAccelerationUnsimplified).cwiseQuotient(currentAccelerationUnsimplified)).transpose()<<std::endl;

//        std::cout<<"mu sun: "<<bodyMap.at( "Sun" )->getGravityFieldModel()->getGravitationalParameter()<<std::endl;
//        std::cout<<"state sun: "<<bodyMap.at( "Sun" )->getState().transpose()<<std::endl;
//        std::cout<<"state mercury: "<<bodyMap.at( "Mercury" )->getState().transpose()<<std::endl;
//        std::cout<<"acceleration: "<<accelerationModel->getAcceleration().transpose()<<std::endl;

        Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );

        Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ) );

        Eigen::MatrixXd partialWrtMercuryPosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtMercuryPosition.block( 0, 0, 3, 3 ) );

        Eigen::MatrixXd partialWrtMercuryVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtMercuryVelocity.block( 0, 0, 3, 3 ) );

        Eigen::Vector3d partialWrtSunGravitationalParameter = accelerationPartial->wrtParameter(
                    gravitationalParameterParameter );

        Eigen::Vector3d partialWrtGamma = accelerationPartial->wrtParameter( ppnParameterGamma );
        Eigen::Vector3d partialWrtBeta = accelerationPartial->wrtParameter( ppnParameterBeta );
        Eigen::Vector3d partialWrtAlpha1 = accelerationPartial->wrtParameter( ppnParameterAlpha1 );
        Eigen::Vector3d partialWrtAlpha2 = accelerationPartial->wrtParameter( ppnParameterAlpha2 );

        Eigen::Vector3d partialWrtNordtvedtParameter = accelerationPartial->wrtParameter( nordtvedtParameter );

        // Declare perturbations in position for numerical partial/
        Eigen::Vector3d positionPerturbation;
        positionPerturbation << 10.0, 10.0, 10.0;
        Eigen::Vector3d velocityPerturbation;
        velocityPerturbation << 1.0, 1.0, 1.0;


        // The partials for position and gravitational parameter have been omitted here
        // In the acceleration and partial calculation long values are used because numerical precision of double's is not sufficient
        // The calculateAccelerationWithStatePartials does not have such functionality
        // To replace the unit test, SEPPartialUnitTest.py is available in python,
        // where it is shown that the partial and central difference matches with tolerance of approximately 10^-15

        //        Eigen::Matrix3d testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
        //                    sunStateSetFunction, accelerationModel, sun->getState( ), positionPerturbation, 0 );
        //        Eigen::Vector3d testPartialWrtSunGravitationalParameter = calculateAccelerationWrtParameterPartials(
        //                    gravitationalParameterParameter, accelerationModel, 1.0E19 );


        Eigen::Matrix3d testPartialWrtSunPosition;
        Eigen::Vector3d testPartialWrtSunGravitationalParameter;
        if (useNordtvedtConstraintFunction()){
            testPartialWrtSunPosition
                    << -5.254917872395E-27, 2.41689255306E-26, 1.45063710477E-26,
                       2.4168925530635E-26, 8.8214542175E-28, -7.3907355489E-27,
                       1.4506371047675E-26, -7.39073554885E-27, 4.37277245065E-27;
            testPartialWrtSunGravitationalParameter
                    << 6.9104839E-50, 9.1491455E-50, 4.1317575E-50;
        } else{
            testPartialWrtSunPosition
                    << -3.9411884042995E-26, 1.812669414795E-25, 1.087977828576E-25,
                       1.81266941479565E-25, 6.6160906632E-27, -5.543051661675E-26,
                       1.087977828576E-25, -5.54305166168E-26, 3.279579337975E-26;
            testPartialWrtSunGravitationalParameter
                    << 3.887147085E-48, 5.14639427E-48, 2.32411352E-48;
        }

        Eigen::Matrix3d testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                    std::bind( &Body::setState, bodyMap.at( "Sun" ), std::placeholders::_1  ),
                    accelerationModel,
                    bodyMap.at( "Sun" )->getState( ),
                    velocityPerturbation, 3 );

        Eigen::Matrix3d testPartialWrtMercuryPosition =
                -1.0 * testPartialWrtSunPosition;

        Eigen::Matrix3d testPartialWrtMercuryVelocity = calculateAccelerationWrtStatePartials(
                    std::bind( &Body::setState, bodyMap.at( "Mercury" ), std::placeholders::_1  ),
                    accelerationModel,
                    bodyMap.at( "Mercury" )->getState( ),
                    velocityPerturbation, 3 );

        Eigen::Vector3d testPartialWrtPpnParameterGamma = calculateAccelerationWrtParameterPartials(
                    ppnParameterGamma, accelerationModel, 100.0 );
        Eigen::Vector3d testPartialWrtPpnParameterBeta = calculateAccelerationWrtParameterPartials(
                    ppnParameterBeta, accelerationModel, 100.0 );

        Eigen::Vector3d testPartialWrtPpnParameterAlpha1 = calculateAccelerationWrtParameterPartials(
                    ppnParameterAlpha1, accelerationModel, 100.0 );
        Eigen::Vector3d testPartialWrtPpnParameterAlpha2 = calculateAccelerationWrtParameterPartials(
                    ppnParameterAlpha2, accelerationModel, 100.0 );

        Eigen::Vector3d testPartialWrtNordtvedtParameter = calculateAccelerationWrtParameterPartials(
                    nordtvedtParameter, accelerationModel, 1.0E-4 );


        // Compare numerical and analytical results.
        double stateTolerance = 1.0E-3; //orginally E-7?
        double parameterTolerance = 1.0E-4; //originally E-8


        std::cout<<"wrt sun position:"
                 <<std::endl<<partialWrtSunPosition
                 <<std::endl<<testPartialWrtSunPosition
                 <<std::endl<<partialWrtSunPosition-testPartialWrtSunPosition
                 <<std::endl;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                           partialWrtSunPosition, stateTolerance);

        std::cout<<"wrt sun velocity:"
                 <<std::endl<<partialWrtSunVelocity
                 <<std::endl<<testPartialWrtSunVelocity
                 <<std::endl<<partialWrtSunVelocity-testPartialWrtSunVelocity
                 <<std::endl;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                           partialWrtSunVelocity, stateTolerance);

        std::cout<<"wrt mercury position:"
                 <<std::endl<<partialWrtMercuryPosition
                 <<std::endl<<testPartialWrtMercuryPosition
                 <<std::endl<<partialWrtMercuryPosition-testPartialWrtMercuryPosition
                 <<std::endl;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMercuryPosition,
                                           partialWrtMercuryPosition, stateTolerance );

        std::cout<<"wrt mercury velocity:"
                 <<std::endl<<partialWrtMercuryVelocity
                 <<std::endl<<testPartialWrtMercuryVelocity
                 <<std::endl<<partialWrtMercuryVelocity-testPartialWrtMercuryVelocity
                 <<std::endl;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMercuryVelocity,
                                           partialWrtMercuryVelocity, stateTolerance );

//        std::cout<<"wrt grav parameter:"
//                 <<std::endl<<partialWrtSunGravitationalParameter.transpose()
//                 <<std::endl<<testPartialWrtSunGravitationalParameter.transpose()
//                 <<std::endl<<partialWrtSunGravitationalParameter.transpose()-testPartialWrtSunGravitationalParameter.transpose()
//                 <<std::endl;

//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunGravitationalParameter,
//                                           partialWrtSunGravitationalParameter, parameterTolerance );

        std::cout<<"wrt gamma:"
                 <<std::endl<<partialWrtGamma.transpose()
                 <<std::endl<<testPartialWrtPpnParameterGamma.transpose()
                 <<std::endl<<partialWrtGamma.transpose()-testPartialWrtPpnParameterGamma.transpose()
                 <<std::endl;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterGamma, partialWrtGamma, parameterTolerance );

        std::cout<<"wrt beta:"
                 <<std::endl<<partialWrtBeta.transpose()
                 <<std::endl<<testPartialWrtPpnParameterBeta.transpose()
                 <<std::endl<<partialWrtBeta.transpose()-testPartialWrtPpnParameterBeta.transpose()
                 <<std::endl;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterBeta, partialWrtBeta, parameterTolerance );

        std::cout<<"wrt alpha1:"
                 <<std::endl<<partialWrtAlpha1.transpose()
                 <<std::endl<<testPartialWrtPpnParameterAlpha1.transpose()
                 <<std::endl<<partialWrtAlpha1.transpose()-testPartialWrtPpnParameterAlpha1.transpose()
                 <<std::endl;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterAlpha1, partialWrtAlpha1, parameterTolerance );

        std::cout<<"wrt alpha2:"
                 <<std::endl<<partialWrtAlpha2.transpose()
                 <<std::endl<<testPartialWrtPpnParameterAlpha2.transpose()
                 <<std::endl<<partialWrtAlpha2.transpose()-testPartialWrtPpnParameterAlpha2.transpose()
                 <<std::endl;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPpnParameterAlpha2, partialWrtAlpha2, parameterTolerance );

        if (useNordtvedtConstraintFunction() == false){
            std::cout<<"wrt nordtvedt parameter:"
                    <<std::endl<<partialWrtNordtvedtParameter.transpose()
                    <<std::endl<<testPartialWrtNordtvedtParameter.transpose()
                    <<std::endl<<partialWrtNordtvedtParameter.transpose()-testPartialWrtNordtvedtParameter.transpose()
                    <<std::endl;

           TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtNordtvedtParameter, partialWrtNordtvedtParameter, parameterTolerance );
        }

    }

}




BOOST_AUTO_TEST_CASE( testTimeVaryingGravitationalParameterPartial )
{
    // Create earth and vehicle bodies.
    std::shared_ptr< Body > sun = std::make_shared< Body >( );
    std::shared_ptr< Body > mercury = std::make_shared< Body >( );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            std::bind( &Body::setState, sun, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > mercuryStateSetFunction =
            std::bind( &Body::setState, mercury, std::placeholders::_1  );
    std::function< Eigen::Vector6d( ) > sunStateGetFunction =
            std::bind( &Body::getState, sun );
    std::function< Eigen::Vector6d( ) > mercuryStateGetFunction =
            std::bind( &Body::getState, mercury);

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    double currentTime = 1.0E8;

    // Set vehicle and earth state.
    sun->setState( getBodyCartesianStateAtEpoch(  "Sun", "SSB", "J2000", "NONE", currentTime ) );
    mercury->setState( getBodyCartesianStateAtEpoch(  "Mercury", "SSB", "J2000", "NONE", currentTime ) );

    NamedBodyMap bodyMap;
    bodyMap[ "Sun" ] = sun;
    bodyMap[ "Mercury" ] = mercury;

    // Create gravity field.
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );
    std::shared_ptr< gravitation::GravityFieldModel > sunGravityField =
            createGravityFieldModel( gravityFieldSettings, "Sun", bodyMap );
    sun->setGravityFieldModel( sunGravityField );
    std::shared_ptr< gravitation::GravityFieldModel > mercuryGravityField =
            createGravityFieldModel( gravityFieldSettings, "Mercury", bodyMap );
    sun->setGravityFieldModel( mercuryGravityField );


    // Create acceleration model.
    ppnParameterSet->setTimeVaryingGravitationalParameter(1.0E-13);

    std::function< double( ) > timeVaryingGravitationalParameterFunction =
            std::bind( &PPNParameterSet::getTimeVaryingGravitationalParameter, ppnParameterSet );

    std::shared_ptr< TimeVaryingGravitationalParameterAcceleration > accelerationModel =
            std::make_shared< TimeVaryingGravitationalParameterAcceleration >
            ( std::bind( &Body::getState, mercury ),
              std::bind( &Body::getState, sun ),
              std::bind( &GravityFieldModel::getGravitationalParameter, sunGravityField ),
              timeVaryingGravitationalParameterFunction );


    // Create acceleration partial object.
    std::shared_ptr< TimeVaryingGravitationalParameterPartial > accelerationPartial =
            std::make_shared< TimeVaryingGravitationalParameterPartial >(
                accelerationModel, "Mercury", "Sun" );

    // Create parameter objects.
    std::shared_ptr< EstimatableParameter< double > > gravitationalParameterParameter = std::make_shared<
            GravitationalParameter >( sunGravityField, "Sun" );
    std::shared_ptr< EstimatableParameter< double > > timeVaryingGravitationalParameter = std::make_shared<
            TimeVaryingGravitationalParameter >( ppnParameterSet );

    // Calculate analytical partials.
    accelerationModel->updateMembers( currentTime );
    accelerationPartial->update( currentTime );

    Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtMercuryPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtMercuryPosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtMercuryVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtMercuryVelocity.block( 0, 0, 3, 3 ) );

    Eigen::Vector3d partialWrtSunGravitationalParameter = accelerationPartial->wrtParameter(
                gravitationalParameterParameter );

    Eigen::Vector3d partialWrtTimeVaryingGravitationalParameter = accelerationPartial->wrtParameter(
                timeVaryingGravitationalParameter );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10.0, 10.0, 10.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0, 1.0;

    // Calculate numerical partials.
    Eigen::Matrix3d testPartialWrtMercuryPosition = calculateAccelerationWrtStatePartials(
                mercuryStateSetFunction, accelerationModel, mercury->getState( ), positionPerturbation, 0, emptyFunction, currentTime );
    Eigen::Matrix3d testPartialWrtMercuryVelocity = calculateAccelerationWrtStatePartials(
                mercuryStateSetFunction, accelerationModel, mercury->getState( ), velocityPerturbation, 3, emptyFunction, currentTime );

    Eigen::Matrix3d testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, accelerationModel, sun->getState( ), positionPerturbation, 0, emptyFunction, currentTime );
    Eigen::Matrix3d testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, accelerationModel, sun->getState( ), velocityPerturbation, 3, emptyFunction, currentTime );

    Eigen::Vector3d testPartialWrtSunGravitationalParameter = calculateAccelerationWrtParameterPartials(
                gravitationalParameterParameter, accelerationModel, 1.0E10, emptyFunction, currentTime );

    Eigen::Vector3d testPartialWrtTimeVaryingGravitationalParameter = calculateAccelerationWrtParameterPartials(
                timeVaryingGravitationalParameter, accelerationModel, 1.0E-14, emptyFunction, currentTime );


    // Compare numerical and analytical results.
    double stateTolerance = 1.0E-3; //orginally E-7?
    double parameterTolerance = 1.0E-4; //originally E-8


    std::cout<<"wrt sun position:"
             <<std::endl<<partialWrtSunPosition
             <<std::endl<<testPartialWrtSunPosition
             <<std::endl<<partialWrtSunPosition-testPartialWrtSunPosition
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                       partialWrtSunPosition, stateTolerance);

    std::cout<<"wrt sun velocity:"
             <<std::endl<<partialWrtSunVelocity
             <<std::endl<<testPartialWrtSunVelocity
             <<std::endl<<partialWrtSunVelocity-testPartialWrtSunVelocity
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                       partialWrtSunVelocity, stateTolerance);

    std::cout<<"wrt mercury position:"
             <<std::endl<<partialWrtMercuryPosition
             <<std::endl<<testPartialWrtMercuryPosition
             <<std::endl<<partialWrtMercuryPosition-testPartialWrtMercuryPosition
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMercuryPosition,
                                       partialWrtMercuryPosition, stateTolerance );

    std::cout<<"wrt mercury velocity:"
             <<std::endl<<partialWrtMercuryVelocity
             <<std::endl<<testPartialWrtMercuryVelocity
             <<std::endl<<partialWrtMercuryVelocity-testPartialWrtMercuryVelocity
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMercuryVelocity,
                                       partialWrtMercuryVelocity, stateTolerance );

    std::cout<<"wrt grav parameter:"
             <<std::endl<<partialWrtSunGravitationalParameter.transpose()
             <<std::endl<<testPartialWrtSunGravitationalParameter.transpose()
             <<std::endl<<partialWrtSunGravitationalParameter.transpose()-testPartialWrtSunGravitationalParameter.transpose()
             <<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunGravitationalParameter,
                                       partialWrtSunGravitationalParameter, parameterTolerance );

    std::cout<<"wrt time varying gravitational parameter:"
             <<std::endl<<partialWrtTimeVaryingGravitationalParameter.transpose()
             <<std::endl<<testPartialWrtTimeVaryingGravitationalParameter.transpose()
             <<std::endl<<partialWrtTimeVaryingGravitationalParameter.transpose()-testPartialWrtTimeVaryingGravitationalParameter.transpose()
             <<std::endl;


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtTimeVaryingGravitationalParameter,
                                       partialWrtTimeVaryingGravitationalParameter, parameterTolerance );


}


BOOST_AUTO_TEST_CASE( testEmpiricalAccelerationPartial )
{

    // Create earth and vehicle bodies.
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    std::shared_ptr< Body > vehicle = std::make_shared< Body >( );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, earth, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            std::bind( &Body::setState, vehicle, std::placeholders::_1  );
    std::function< Eigen::Vector6d( ) > earthStateGetFunction =
            std::bind( &Body::getState, earth );
    std::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
            std::bind( &Body::getState, vehicle );

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Set vehicle and earth state.
    earth->setState( getBodyCartesianStateAtEpoch(  "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );
    Eigen::Vector6d vehicleKeplerElements;
    vehicleKeplerElements << 6378.0E3 + 249E3, 0.0004318, convertDegreesToRadians( 96.5975 ),
            convertDegreesToRadians( 217.6968 ), convertDegreesToRadians( 268.2663 ), convertDegreesToRadians( 142.3958 );
    vehicle->setState( earth->getState( ) + convertKeplerianToCartesianElements( vehicleKeplerElements,
                                                                                 getBodyGravitationalParameter( "Earth" ) ) );


    NamedBodyMap bodyMap;
    bodyMap[ "Vehicle" ] = vehicle;
    bodyMap[ "Earth" ] = earth;

    // Create gravity field.
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );
    std::shared_ptr< gravitation::GravityFieldModel > earthGravityField =
            createGravityFieldModel( gravityFieldSettings, "Earth", bodyMap );
    earth->setGravityFieldModel( earthGravityField );

    // Create rotation model
    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > simpleRotationalEphemeris =
            std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000" , "IAU_Earth", 0.0 ),
                2.0 * mathematical_constants::PI / 86400.0,
                1.0E7,
                "ECLIPJ2000" , "IAU_Earth" );
    earth->setRotationalEphemeris( simpleRotationalEphemeris );
    earth->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );

    // Define empirical acceleration
    Eigen::Vector3d constantAcceleration = 0.0 * Eigen::Vector3d( 0.038, -0.7528, 0.00752 );
    Eigen::Vector3d sineAcceleration = Eigen::Vector3d( 0.984, 0.0427, -0.0764238 );
    Eigen::Vector3d cosineAcceleration = Eigen::Vector3d( -0.0024785, 1.839, -0.73288 );

    // Create acceleration model.
    std::shared_ptr< EmpiricalAcceleration > accelerationModel =
            std::make_shared< EmpiricalAcceleration >
            ( constantAcceleration, sineAcceleration, cosineAcceleration,
              std::bind( &Body::getState, vehicle ),
              std::bind( &GravityFieldModel::getGravitationalParameter, earthGravityField ),
              std::bind( &Body::getState, earth ) );

    // Create acceleration partial object.
    std::shared_ptr< EmpiricalAccelerationPartial > accelerationPartial = std::make_shared< EmpiricalAccelerationPartial >(
                accelerationModel, "Vehicle", "Earth" );

    // Define list of empirical accelerations w.r.t. which partials are to be computed
    std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > allEmpiricalShapesVector;
    allEmpiricalShapesVector.push_back( basic_astrodynamics::constant_empirical );
    allEmpiricalShapesVector.push_back( basic_astrodynamics::cosine_empirical );

    std::map< basic_astrodynamics::EmpiricalAccelerationComponents, std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > >
            empiricalComponentsToEstimate;
    empiricalComponentsToEstimate[ basic_astrodynamics::radial_empirical_acceleration_component ] = allEmpiricalShapesVector;
    empiricalComponentsToEstimate[ basic_astrodynamics::along_track_empirical_acceleration_component ] = allEmpiricalShapesVector;

    allEmpiricalShapesVector.push_back( basic_astrodynamics::sine_empirical );
    empiricalComponentsToEstimate[ basic_astrodynamics::across_track_empirical_acceleration_component ] = allEmpiricalShapesVector;

    // Create time-independent empirical acceleration object.
    std::shared_ptr< EmpiricalAccelerationCoefficientsParameter > empiricalAccelerationParameter = std::make_shared<
            EmpiricalAccelerationCoefficientsParameter >(
    std::vector< std::shared_ptr< EmpiricalAcceleration > >( { accelerationModel } ), "Vehicle", "Earth",
                empiricalComponentsToEstimate );

    {
        // Calculate analytical partials.
        accelerationModel->updateMembers( );
        accelerationPartial->update( );
        Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtEmpiricalCoefficients;
        accelerationPartial->wrtEmpiricalAccelerationCoefficient(
                    empiricalAccelerationParameter, partialWrtEmpiricalCoefficients );

        // Declare perturbations in position for numerical partial
        Eigen::Vector3d positionPerturbation;
        positionPerturbation << 1.0, 1.0, 1.0;
        Eigen::Vector3d velocityPerturbation;
        velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;
        int parameterSize = empiricalAccelerationParameter->getParameterSize( );
        Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Constant( parameterSize, 1.0E-5 );

        // Calculate numerical partials.
        Eigen::MatrixXd testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, accelerationModel, vehicle->getState( ), 10.0 * positionPerturbation, 0 );
        Eigen::MatrixXd testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, accelerationModel, vehicle->getState( ), velocityPerturbation, 3 );
        Eigen::MatrixXd testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                    earthStateSetFunction, accelerationModel, earth->getState( ), positionPerturbation, 0 );
        Eigen::MatrixXd testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                    earthStateSetFunction, accelerationModel, earth->getState( ), velocityPerturbation, 3 );
        Eigen::MatrixXd testPartialWrtEmpiricalCoefficients = calculateAccelerationWrtParameterPartials(
                    empiricalAccelerationParameter, accelerationModel, parameterPerturbation );

        // Compare numerical and analytical results.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                           partialWrtEarthPosition, 1.0e-3 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                           partialWrtEarthVelocity, 1.0e-3 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                           partialWrtVehiclePosition, 1.0e-3 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                           partialWrtVehicleVelocity, 1.0e-3 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEmpiricalCoefficients,
                                           partialWrtEmpiricalCoefficients, 1.0e-3 );
    }

    // Define arc split times of arc-wise empirical accelerations
    std::vector< double > arcStartTimes;
    arcStartTimes.push_back( 0.0 );
    arcStartTimes.push_back( 1.0E4 );
    arcStartTimes.push_back( 2.0E4 );
    arcStartTimes.push_back( 3.0E4 );
    arcStartTimes.push_back( 7.0E4 );
    std::shared_ptr< ArcWiseEmpiricalAccelerationCoefficientsParameter > arcWiseEmpiricalAccelerationParameter = std::make_shared<
            ArcWiseEmpiricalAccelerationCoefficientsParameter >(
                std::vector< std::shared_ptr< EmpiricalAcceleration > >( { accelerationModel } ), "Vehicle", "Earth",
                empiricalComponentsToEstimate, arcStartTimes );

    // Define list of times at which to test empirical acceleration
    std::vector< double > evaluationTimes;

    evaluationTimes.push_back( 1.0 );
    evaluationTimes.push_back( 0.5E4 );
    evaluationTimes.push_back( 1.0E4 - 1.0 );
    evaluationTimes.push_back( 1.0E4 + 1.0 );
    evaluationTimes.push_back( 1.2E4 );
    evaluationTimes.push_back( 3.5E4 );
    evaluationTimes.push_back( 1.0E5 );

    Eigen::VectorXd accelerationPerturbationVector = Eigen::VectorXd::Zero( arcWiseEmpiricalAccelerationParameter->getParameterSize( ) );
    for( int i = 0; i < arcWiseEmpiricalAccelerationParameter->getParameterSize( ); i++  )
    {
        accelerationPerturbationVector( i ) = 1.0E-6;
    }

    // Iterate over all test times.
    for( unsigned int i = 0; i < evaluationTimes.size( ); i++ )
    {
        // Update models to current time
        accelerationModel->resetTime( TUDAT_NAN );
        accelerationModel->updateMembers( evaluationTimes.at( i ) );
        accelerationPartial->update( evaluationTimes.at( i ) );

        // Compute analytical partials
        Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtEmpiricalCoefficients = accelerationPartial->wrtParameter(
                    arcWiseEmpiricalAccelerationParameter );

        // Set numerical partial settings
        Eigen::Vector3d positionPerturbation;
        positionPerturbation << 1.0, 1.0, 1.0;
        Eigen::Vector3d velocityPerturbation;
        velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;
        int parameterSize = empiricalAccelerationParameter->getParameterSize( );
        Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Constant( parameterSize, 1.0E-5 );

        // Calculate numerical partials.
        Eigen::MatrixXd testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, accelerationModel, vehicle->getState( ), positionPerturbation, 0, emptyFunction,
                    evaluationTimes.at( i ) );
        Eigen::MatrixXd testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                    vehicleStateSetFunction, accelerationModel, vehicle->getState( ), velocityPerturbation, 3, emptyFunction,
                    evaluationTimes.at( i ) );
        Eigen::MatrixXd testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                    earthStateSetFunction, accelerationModel, earth->getState( ), positionPerturbation, 0, emptyFunction,
                    evaluationTimes.at( i ) );
        Eigen::MatrixXd testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                    earthStateSetFunction, accelerationModel, earth->getState( ), velocityPerturbation, 3, emptyFunction,
                    evaluationTimes.at( i ) );
        Eigen::MatrixXd testPartialWrtEmpiricalCoefficients = calculateAccelerationWrtParameterPartials(
                    arcWiseEmpiricalAccelerationParameter, accelerationModel, accelerationPerturbationVector,
                    emptyFunction, evaluationTimes.at( i ) );


        //Compare numerical and analytical results.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition,
                                           partialWrtEarthPosition, 1.0e-5 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity,
                                           partialWrtEarthVelocity, 1.0e-6 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                           partialWrtVehiclePosition, 1.0e-5 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                           partialWrtVehicleVelocity, 1.0e-6 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEmpiricalCoefficients,
                                           partialWrtEmpiricalCoefficients, 1.0e-6 );
    }
}

BOOST_AUTO_TEST_CASE( testDirectDissipationAccelerationPartial )
{

    // Create bodies
    std::shared_ptr< Body > jupiter = std::make_shared< Body >( );
    std::shared_ptr< Body > io = std::make_shared< Body >( );

    // Create links to set and get state functions of bodies.
    std::function< void( Eigen::Vector6d ) > ioStateSetFunction =
            std::bind( &Body::setState, io, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > jupiterStateSetFunction =
            std::bind( &Body::setState, jupiter, std::placeholders::_1  );
    std::function< Eigen::Vector6d( ) > ioStateGetFunction =
            std::bind( &Body::getState, io );
    std::function< Eigen::Vector6d( ) > jupiterStateGetFunction =
            std::bind( &Body::getState, jupiter );

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Set state.
    jupiter->setState( Eigen::Vector6d::Zero( ) );
    Eigen::Vector6d ioKeplerElements =
            ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.001, 2.0, 3.0, 0.4 ).finished( );
    io->setState( convertKeplerianToCartesianElements(
                      ioKeplerElements,
                      getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ) ) );


    NamedBodyMap bodyMap;
    bodyMap[ "Io" ] = io;
    bodyMap[ "Jupiter" ] = jupiter;


    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );

    // Create jupiter gravity field.
    std::shared_ptr< GravityFieldSettings > jupiterGravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >
            ( getBodyGravitationalParameter( "Jupiter" ), getAverageRadius( "Jupiter" ),
              cosineCoefficients, sineCoefficients, "IAU_Jupiter" );
    std::shared_ptr< gravitation::GravityFieldModel > jupiterGravityField =
            createGravityFieldModel( jupiterGravityFieldSettings, "Jupiter", bodyMap );
    jupiter->setGravityFieldModel( jupiterGravityField );

    // Create io gravity field.
    std::shared_ptr< GravityFieldSettings > ioGravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >
            ( getBodyGravitationalParameter( "Io" ), getAverageRadius( "Io" ),
              cosineCoefficients, sineCoefficients, "IAU_Io" );
    std::shared_ptr< gravitation::GravityFieldModel > ioGravityField =
            createGravityFieldModel( ioGravityFieldSettings, "Io", bodyMap );
    io->setGravityFieldModel( ioGravityField );

    // Create rotation model
    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > simpleRotationalEphemeris =
            std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ),
                2.0 * mathematical_constants::PI / ( 9.925 * 3600.0 ), 0.0,
                "ECLIPJ2000" , "IAU_Jupiter" );
    jupiter->setRotationalEphemeris( simpleRotationalEphemeris );
    jupiter->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );

    double loveNumber = 0.1;
    double timeLag = 100.0;
    for( unsigned int useRadialTerm = 0; useRadialTerm < 2; useRadialTerm++ )
    {
        for( unsigned int usePlanetTide = 0; usePlanetTide < 2; usePlanetTide++ )
        {

            // Create acceleration model.
            std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > accelerationModel =
                    std::dynamic_pointer_cast<  gravitation::DirectTidalDissipationAcceleration >(
                        simulation_setup::createAccelerationModel(
                            io, jupiter, std::make_shared< simulation_setup::DirectTidalDissipationAccelerationSettings >(
                                loveNumber, timeLag, useRadialTerm, usePlanetTide ) , "Io", "Jupiter" ) );

            // Create acceleration partial object.
            std::shared_ptr< acceleration_partials::DirectTidalDissipationAccelerationPartial > accelerationPartial =
                    std::make_shared< acceleration_partials::DirectTidalDissipationAccelerationPartial >(
                        accelerationModel, "Io", "Jupiter" );

            // Create gravitational parameter object.
            std::shared_ptr< EstimatableParameter< double > > jupiterGravitationalParameterParameter = std::make_shared<
                    GravitationalParameter >( jupiterGravityField, "Jupiter" );
            std::shared_ptr< EstimatableParameter< double > > ioGravitationalParameterParameter = std::make_shared<
                    GravitationalParameter >( ioGravityField, "Io" );
            std::shared_ptr< EstimatableParameter< double > > tidalTimeLagParameter = std::make_shared< DirectTidalTimeLag >(
                        std::vector< std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > >{ accelerationModel },
                        usePlanetTide ? "Jupiter" : "Io" );


            {
                // Calculate analytical partials.
                accelerationModel->updateMembers( );
                accelerationPartial->update( );
                Eigen::MatrixXd partialWrtJupiterPosition = Eigen::Matrix3d::Zero( );
                accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtJupiterPosition.block( 0, 0, 3, 3 ) );
                Eigen::MatrixXd partialWrtJupiterVelocity = Eigen::Matrix3d::Zero( );
                accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtJupiterVelocity.block( 0, 0, 3, 3 ) );
                Eigen::MatrixXd partialWrtIoPosition = Eigen::Matrix3d::Zero( );
                accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtIoPosition.block( 0, 0, 3, 3 ) );
                Eigen::MatrixXd partialWrtIoVelocity = Eigen::Matrix3d::Zero( );
                accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtIoVelocity.block( 0, 0, 3, 3 ) );
                Eigen::MatrixXd partialWrtJupiterGravitationalParameter = accelerationPartial->wrtParameter(
                            jupiterGravitationalParameterParameter );
                Eigen::MatrixXd partialWrtIoGravitationalParameter = accelerationPartial->wrtParameter(
                            ioGravitationalParameterParameter );
                Eigen::MatrixXd partialWrtTidalTimeLag = accelerationPartial->wrtParameter(
                            tidalTimeLagParameter );

                // Declare perturbations in position for numerical partial
                Eigen::Vector3d positionPerturbation;
                positionPerturbation << 10.0, 10.0, 10.0;
                Eigen::Vector3d velocityPerturbation;
                velocityPerturbation << 1.0E-1, 1.0E-1, 1.0E-1;
                double jupiterGravityFieldPerturbation = 1.0E8;
                double ioGravityFieldPerturbation = 1.0E8;
                double timelagPerturbation = 1.0;;

                // Calculate numerical partials.
                Eigen::MatrixXd testPartialWrtIoPosition = calculateAccelerationWrtStatePartials(
                            ioStateSetFunction, accelerationModel, io->getState( ), positionPerturbation, 0 );
                Eigen::MatrixXd testPartialWrtIoVelocity = calculateAccelerationWrtStatePartials(
                            ioStateSetFunction, accelerationModel, io->getState( ), velocityPerturbation, 3 );
                Eigen::MatrixXd testPartialWrtJupiterPosition = calculateAccelerationWrtStatePartials(
                            jupiterStateSetFunction, accelerationModel, jupiter->getState( ), positionPerturbation, 0 );
                Eigen::MatrixXd testPartialWrtJupiterVelocity = calculateAccelerationWrtStatePartials(
                            jupiterStateSetFunction, accelerationModel, jupiter->getState( ), velocityPerturbation, 3 );
                Eigen::MatrixXd testPartialWrtJupiterGravitationalParameter = calculateAccelerationWrtParameterPartials(
                            jupiterGravitationalParameterParameter, accelerationModel, jupiterGravityFieldPerturbation );
                Eigen::MatrixXd testPartialWrtIoGravitationalParameter = calculateAccelerationWrtParameterPartials(
                            ioGravitationalParameterParameter, accelerationModel, ioGravityFieldPerturbation );
                Eigen::MatrixXd testPartialWrtTidalTimeLag = calculateAccelerationWrtParameterPartials(
                            tidalTimeLagParameter, accelerationModel, timelagPerturbation );

                // Compare numerical and analytical results.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtJupiterPosition,
                                                   partialWrtJupiterPosition, 1.0e-5 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtJupiterVelocity,
                                                   partialWrtJupiterVelocity, 1.0e-5 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtIoPosition,
                                                   partialWrtIoPosition, 1.0e-5 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtIoVelocity,
                                                   partialWrtIoVelocity, 1.0e-5 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtJupiterGravitationalParameter,
                                                   partialWrtJupiterGravitationalParameter, 1.0e-6 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtIoGravitationalParameter,
                                                   partialWrtIoGravitationalParameter, 1.0e-6 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtTidalTimeLag,
                                                   partialWrtTidalTimeLag, 1.0e-6 );

            }
        }
    }
}


BOOST_AUTO_TEST_CASE( testPanelledRadiationPressureAccelerationPartials )
{
    // Create empty bodies, earth and sun.
    std::shared_ptr< Body > vehicle = std::make_shared< Body >( );
    double vehicleMass = 400.0;
    vehicle->setConstantBodyMass( vehicleMass );
    std::shared_ptr< Body > sun = std::make_shared< Body >( );
    NamedBodyMap bodyMap;
    bodyMap[ "Vehicle" ] = vehicle;
    bodyMap[ "Sun" ] = sun;

    // Load spice kernels.
    tudat::spice_interface::loadStandardSpiceKernels( );

    // Set current state of sun and earth.
    sun->setState( Eigen::Vector6d::Zero( ) );//getBodyCartesianStateAtEpoch( "Sun", "SSB", "J2000", "NONE", 1.0E6 ) );
    vehicle->setState(
                ( Eigen::Vector6d( ) << 1.4E11, 1.0E11, 1.1E11, 0.0, 0.0, 0.0 ).finished( ) );//getBodyCartesianStateAtEpoch(  "Earth", "SSB", "J2000", "NONE", 1.0E6 ) );
    vehicle->setRotationalEphemeris(
                std::make_shared< tudat::ephemerides::SimpleRotationalEphemeris >( 0.2, 0.4, -0.2, 1.0E-5, 0.0, "ECLIPJ2000", "VehicleFixed" ) );
    vehicle->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );

    // Create links to set and get state functions of bodies.
    boost::function< void( Eigen::Vector6d ) > sunStateSetFunction =
            boost::bind( &Body::setState, sun, _1  );
    boost::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            boost::bind( &Body::setState, vehicle, _1  );
    boost::function< Eigen::Vector6d( ) > sunStateGetFunction =
            boost::bind( &Body::getState, sun );
    boost::function< Eigen::Vector6d( ) > vehicleStateGetFunction =
            boost::bind( &Body::getState, vehicle );

    // Create radiation pressure properties of vehicle
    std::vector< double > areas;
    areas.push_back( 1.0 );
    areas.push_back( 3.254 );
    areas.push_back( 8.654 );
    areas.push_back( 1.346 );
    areas.push_back( 10.4783 );
    areas.push_back( 6.4235 );
    areas.push_back( 12.483 );

    std::vector< double > emissivities;
    emissivities.push_back( 0.5 );
    emissivities.push_back( 0.9 );
    emissivities.push_back( 0.6 );
    emissivities.push_back( 0.3 );
    emissivities.push_back( 0.4 );
    emissivities.push_back( 0.2 );
    emissivities.push_back( 0.5 );

    std::vector< double > diffuseReflectionCoefficients;
    diffuseReflectionCoefficients.push_back( 0.1 );
    diffuseReflectionCoefficients.push_back( 0.2 );
    diffuseReflectionCoefficients.push_back( 0.15 );
    diffuseReflectionCoefficients.push_back( 0.02 );
    diffuseReflectionCoefficients.push_back( 0.05 );
    diffuseReflectionCoefficients.push_back( 0.12 );
    diffuseReflectionCoefficients.push_back( 0.08 );

    std::vector< Eigen::Vector3d > panelSurfaceNormals;
    panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitX( ) );
    panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitY( ) );
    panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitZ( ) );
    panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
    panelSurfaceNormals.push_back( Eigen::Vector3d::UnitY( ) );
    panelSurfaceNormals.push_back( Eigen::Vector3d::UnitZ( ) );
    panelSurfaceNormals.push_back( ( Eigen::Vector3d( ) << -0.5 * sqrt( 2.0 ), -0.5 * sqrt( 2.0 ), 0.0 ).finished( ) );


    std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface =
            std::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                createRadiationPressureInterface( std::make_shared< PanelledRadiationPressureInterfaceSettings >(
                        "Sun", areas, emissivities, diffuseReflectionCoefficients, panelSurfaceNormals ), "Vehicle", bodyMap ) );

    radiationPressureInterface->updateInterface( 0.0 );
    vehicle->setRadiationPressureInterface( "Sun", radiationPressureInterface );

    // Create acceleration model.
    std::shared_ptr< PanelledRadiationPressureAcceleration > accelerationModel =
            std::make_shared< PanelledRadiationPressureAcceleration >(
                radiationPressureInterface, std::bind( &Body::getBodyMass, vehicle ) );
    accelerationModel->updateMembers( );

    // Create partial-calculating object.
    std::shared_ptr< PanelledRadiationPressurePartial > accelerationPartial =
            std::make_shared< PanelledRadiationPressurePartial >
            ( accelerationModel, radiationPressureInterface, "Vehicle", "Sun" );

//    std::vector< int > panelIndices1 = boost::assign::list_of( 0 )( 6 );
//    std::vector< int > panelIndices2 = boost::assign::list_of( 2 );

    //    std::vector< std::vector< int > > panelIndices;
    //    panelIndices.push_back( panelIndices2 );
    //    panelIndices.push_back( panelIndices1 );

    //    boost::shared_ptr< EstimatableParameterSettings > panelEmissivitiesSettings =
    //            boost::make_shared< PanelRadiationEmissivitiesParameterSettings >( "Vehicle", panelIndices );
    //    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterSettingsVector;
    //    parameterSettingsVector.push_back( panelEmissivitiesSettings );
    //    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > estimatableParameters = createParametersToEstimate(
    //                parameterSettingsVector, bodyMap );
    //    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > panelEmissivitiesParameter =
    //            estimatableParameters->getVectorParameters( ).begin( )->second;

    // Calculate analytical partials.
    accelerationPartial->update( );
    Eigen::MatrixXd partialWrtSunPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtSunPosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtSunVelocity.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ) );

    //Eigen::MatrixXd partialWrtEmissivities = accelerationPartial->wrtParameter( panelEmissivitiesParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtSunVelocity = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd testPartialWrtEmissivities = Eigen::MatrixXd::Zero( 3, 2 );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation<< 1000.0, 1000.0, 1000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation<< 0.1, 0.1, 0.1;
    //    Eigen::VectorXd emissivityPerturbations = Eigen::VectorXd::Zero( 2 );
    //    emissivityPerturbations( 0 ) = 1.0;
    //    emissivityPerturbations( 1 ) = 1.0;

    // Calculate numerical partials.
    std::function< void( ) > updateFunction = std::bind( &RadiationPressureInterface::updateInterface,
                                                             radiationPressureInterface, 0.0 );
    testPartialWrtSunPosition = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, accelerationModel, sun->getState( ), positionPerturbation, 0, updateFunction );
    testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, vehicle->getState( ), positionPerturbation, 0, updateFunction );
    testPartialWrtSunVelocity = calculateAccelerationWrtStatePartials(
                sunStateSetFunction, accelerationModel, sun->getState( ),velocityPerturbation, 3, updateFunction );
    testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, accelerationModel, vehicle->getState( ), velocityPerturbation, 3, updateFunction );
    //    testPartialWrtEmissivities = calculateAccelerationWrtParameterPartials(
    //                panelEmissivitiesParameter, accelerationModel, emissivityPerturbations );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunPosition,
                                       partialWrtSunPosition, 1.0e-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSunVelocity,
                                       partialWrtSunVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition,
                                       partialWrtVehiclePosition, 1.0e-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity,
                                       partialWrtVehicleVelocity, std::numeric_limits< double >::epsilon( ) );
    //    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEmissivities,
    //                                       partialWrtEmissivities, 1.0e-14 );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




