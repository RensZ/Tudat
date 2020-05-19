/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/relativisticAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function to compute partial of Schwarzschild acceleration correction w.r.t. position of body undergoing acceleration
void computePartialOfSchwarzschildAccelerationCorrectionWrtPosition(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter,
        const double ppnParameterGamma,
        const double ppnParameterBeta )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = 2.0 * ( ppnParameterGamma + ppnParameterBeta ) * gravitationalParameter / distance * (
                Eigen::Matrix3d::Identity( ) - position * position.transpose( ) / ( distance * distance ) );
    partialMatrix -= ppnParameterGamma * velocity.dot( velocity ) * Eigen::Matrix3d::Identity( );
    partialMatrix += 2.0 * ( 1.0 + ppnParameterGamma ) * velocity * velocity.transpose( );
    partialMatrix *= gravitationalParameter *  physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance );
    partialMatrix -= 3.0 * currentAcceleration * position.transpose( ) / ( distance * distance );
}

//! Function to compute partial of Schwarzschild acceleration correction w.r.t. velocity of body undergoing acceleration
void computePartialOfSchwarzschildAccelerationCorrectionWrtVelocity(
        const Eigen::Vector6d& relativeState,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter,
        const double ppnParameterGamma )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = gravitationalParameter * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance ) *
            ( - 2.0 * ppnParameterGamma * position * velocity.transpose( ) +
              2.0 * ( 1.0 + ppnParameterGamma ) * (
                  position.dot( velocity ) * Eigen::Matrix3d::Identity( ) + velocity * position.transpose( ) ) );

}

//! Function to compute partial derivative of Schwarzschild acceleration correction w.r.t. central body gravitational patameter.
void computePartialOfSchwarzschildAccelerationCorrectionWrtGravitationalParameter(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partialMatrix,
        const double ppnParameterGamma,
        const double ppnParameterBeta )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance ) *
            ( -ppnParameterGamma * ( velocity.dot( velocity ) ) * position +
              2.0 * ( 1.0 + ppnParameterGamma ) *
              ( position.dot( velocity ) ) * velocity
              + 4.0 * gravitationalParameter * ( ppnParameterGamma + ppnParameterBeta ) *
              position / distance );
}

//! Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter gamma
void computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterGamma(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partialMatrix )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * gravitationalParameter / ( distance * distance * distance ) * (
                ( 2.0 * gravitationalParameter / distance - velocity.dot( velocity ) ) * position + 2.0 * position.dot( velocity ) * velocity );
    //std::cout<<"g:"<<partialMatrix.transpose()<<std::endl;
}

//! Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter beta
void computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterBeta(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partialMatrix )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    double distance = position.norm( );

    partialMatrix = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * gravitationalParameter / ( distance * distance * distance ) * (
                ( 2.0 * gravitationalParameter / distance ) * position );
    //std::cout<<"b:"<<partialMatrix.transpose()<<std::endl;
}





void computePartialOfSchwarzschildAlphaTermsWrtPosition(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameterBodyExertingAcceleration,
        const double gravitationalParameterBodyUndergoingAcceleration,
        const double ppnParameterAlpha1,
        const double ppnParameterAlpha2)
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );
    double gravitationalParameterRatio =
            (gravitationalParameterBodyUndergoingAcceleration / gravitationalParameterBodyExertingAcceleration);

//    partialMatrix = 2.0 * ( ppnParameterGamma + ppnParameterBeta ) * gravitationalParameter / distance * (
//                Eigen::Matrix3d::Identity( ) - position * position.transpose( ) / ( distance * distance ) );

    partialMatrix = (2.0 + ppnParameterAlpha1)
            * (gravitationalParameterBodyUndergoingAcceleration / distance)
            * (Eigen::Matrix3d::Identity( ) - position * position.transpose() / (distance * distance));
    partialMatrix -= 0.5 * (6.0 + ppnParameterAlpha1 + ppnParameterAlpha2)
            * gravitationalParameterRatio
            * (velocity.dot(velocity))
            * Eigen::Matrix3d::Identity( );

//    partialMatrix += 1.5 * (1.0 + ppnParameterAlpha2)
//            * gravitationalParameterRatio
//            / (distance * distance)
//            * velocity.dot(position) * velocity.dot(position)
//            * ( (3.0 * Eigen::Matrix3d::Identity( ) - 2.0 * position * position.transpose() / (distance * distance)));
    partialMatrix += 1.5 * (1.0 + ppnParameterAlpha2)
            * gravitationalParameterRatio
            * velocity.dot(position)
            / (distance * distance )
            * ( 2.0 * position * velocity.transpose()
                - 2.0 * velocity.dot(position) * position * position.transpose() / (distance * distance)
                + velocity.dot(position) * Eigen::Matrix3d::Identity()
                );

    partialMatrix -= gravitationalParameterRatio
            * (2.0 - ppnParameterAlpha1 + ppnParameterAlpha2)
            * velocity * velocity.transpose();
    partialMatrix *= gravitationalParameterBodyExertingAcceleration
            * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT
            / ( distance * distance * distance );
    partialMatrix -= 3.0 * currentAcceleration * position.transpose( )
            / ( distance * distance );
//    std::cout<<" pos partial: "<<std::endl<<partialMatrix<<std::endl;

}

void computePartialOfSchwarzschildAlphaTermsWrtVelocity(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameterBodyExertingAcceleration,
        const double gravitationalParameterBodyUndergoingAcceleration,
        const double ppnParameterAlpha1,
        const double ppnParameterAlpha2)
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = -1.0 * (6.0 + ppnParameterAlpha1 + ppnParameterAlpha2)
            * position * velocity.transpose( );

    partialMatrix += 3.0 * (1.0 + ppnParameterAlpha2)
            * velocity.dot(position) * position * position.transpose()
            / (distance * distance);

    partialMatrix -= (2.0 - ppnParameterAlpha1 + ppnParameterAlpha2)
            * ( position.dot(velocity) * Eigen::Matrix3d::Identity( ) + velocity * position.transpose() );

    partialMatrix *= gravitationalParameterBodyUndergoingAcceleration
            * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT
            / ( distance * distance * distance );
//    std::cout<<" vel partial: "<<std::endl<<partialMatrix<<std::endl;
}

void computePartialOfSchwarzschildAlphaTermsWrtGravitationalParameter(
        Eigen::Vector3d& currentAcceleration,
        Eigen::MatrixXd& partialMatrix,
        const double gravitationalParameterBodyExertingAcceleration)
{
    partialMatrix = currentAcceleration / gravitationalParameterBodyExertingAcceleration;
}


void computePartialOfSchwarzschildAlphaTermsWrtPpnParameterAlpha1(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        Eigen::MatrixXd& partialMatrix,
        const double gravitationalParameterBodyExertingAcceleration,
        const double gravitationalParameterBodyUndergoingAcceleration)
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );
    double gravitationalParameterRatio =
            (gravitationalParameterBodyUndergoingAcceleration / gravitationalParameterBodyExertingAcceleration);

    partialMatrix = (
                (gravitationalParameterBodyUndergoingAcceleration / distance )
                - 0.5 * gravitationalParameterRatio * (velocity.dot(velocity))
            ) * position;
    partialMatrix += gravitationalParameterRatio
            * position.dot(velocity)
            * velocity;
    partialMatrix *= gravitationalParameterBodyExertingAcceleration
            * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT
            / ( distance * distance * distance );

    //std::cout<<"a1:"<<partialMatrix.transpose()<<std::endl;
}


void computePartialOfSchwarzschildAlphaTermsWrtPpnParameterAlpha2(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        Eigen::MatrixXd& partialMatrix,
        const double gravitationalParameterBodyUndergoingAcceleration)
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = (
                -0.5 * (velocity.dot(velocity))
                + 1.5 * (velocity.dot(position/ distance)) * (velocity.dot(position/ distance))
            ) * position;
    partialMatrix -= position.dot(velocity) * velocity;
    partialMatrix *= gravitationalParameterBodyUndergoingAcceleration *
            physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT /
            ( distance * distance * distance );

    //std::cout<<"a2:"<<partialMatrix.transpose()<<std::endl;
}





//! Function to compute partial of LenseThirring acceleration correction w.r.t. position of body undergoing acceleration
void computePartialOfLenseThirringAccelerationCorrectionWrtPosition(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter,
        const double ppnParameterGamma)
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    Eigen::Matrix3d velocityCrossMatrix;
    velocityCrossMatrix << 0.0, -velocity(2), velocity(1),
                           velocity(2), 0.0, -velocity(0),
                           -velocity(1), velocity(0), 0.0;

    partialMatrix = position.cross(velocity)
            * position.dot(centralBodyAngularMomentum)
            * -2.0 * position.transpose() / (distance * distance);
    partialMatrix += position.dot(centralBodyAngularMomentum)
            * velocityCrossMatrix.transpose();
    partialMatrix += position.cross(velocity)
            * centralBodyAngularMomentum.transpose();

    partialMatrix *= 3.0 / (distance * distance);

    partialMatrix *= (1.0 + ppnParameterGamma) * gravitationalParameter
            *  physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance );
    partialMatrix -= 3.0 * currentAcceleration * position.transpose( ) / ( distance * distance );

//    partialMatrix = Eigen::Matrix3d::Zero();
}



//! Function to compute partial of LenseThirring acceleration correction w.r.t. velocity of body undergoing acceleration
void computePartialOfLenseThirringAccelerationCorrectionWrtVelocity(
        const Eigen::Vector6d& relativeState,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter,
        const double ppnParameterGamma )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    double distance = position.norm( );

    Eigen::Matrix3d positionCrossMatrix;
    positionCrossMatrix << 0.0, -position(2), position(1),
                           position(2), 0.0, -position(0),
                           -position(1), position(0), 0.0;

    Eigen::Matrix3d angularMomentumCrossMatrix;
    angularMomentumCrossMatrix << 0.0, -centralBodyAngularMomentum(2), centralBodyAngularMomentum(1),
                                  centralBodyAngularMomentum(2), 0.0, -centralBodyAngularMomentum(0),
                                  -centralBodyAngularMomentum(1), centralBodyAngularMomentum(0), 0.0;

    partialMatrix = 3.0
            * position.dot(centralBodyAngularMomentum)
            * positionCrossMatrix
            / (distance * distance);

    partialMatrix += angularMomentumCrossMatrix.transpose();

    partialMatrix *= (1.0 + ppnParameterGamma) * gravitationalParameter
            *  physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance );

//    partialMatrix = Eigen::Matrix3d::Zero();
}

//! Function to compute partial derivative of LenseThirring acceleration correction w.r.t. central body gravitational patameter.
void computePartialOfLenseThirringAccelerationCorrectionWrtGravitationalParameter(
        const Eigen::Vector6d& relativeState,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        Eigen::MatrixXd& partialMatrix,
        const double ppnParameterGamma)
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    double commonCorrectionTerm = (1.0 + ppnParameterGamma)/
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT *
              distance * distance * distance );
    Eigen::Vector3d termInBrackets =
            3.0 / ( distance * distance ) *
            position.cross( velocity ) *
            ( position.dot( centralBodyAngularMomentum ) ) +
            velocity.cross( centralBodyAngularMomentum );

    partialMatrix = commonCorrectionTerm * termInBrackets;
//    partialMatrix = Eigen::Vector3d::Zero();
}



//! Function to compute the partial derivative of LenseThirring acceleration correction w.r.t. PPN parameter gamma
void computePartialOfLenseThirringAccelerationCorrectionWrtPpnParameterGamma(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        Eigen::MatrixXd& partialMatrix )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    double commonCorrectionTerm = gravitationalParameter/
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT *
              distance * distance * distance );
    Eigen::Vector3d termInBrackets =
            3.0 / ( distance * distance ) *
            position.cross( velocity ) *
            ( position.dot( centralBodyAngularMomentum ) ) +
            velocity.cross( centralBodyAngularMomentum );

    partialMatrix = commonCorrectionTerm * termInBrackets;
//    partialMatrix = Eigen::Vector3d::Zero();
}




//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
RelativisticAccelerationPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    // Create partial function if parameter is central body gravitational parameter
    if( parameter->getParameterName( ).second.first == acceleratingBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case estimatable_parameters::gravitational_parameter:
            partialFunction = std::bind( &RelativisticAccelerationPartial::wrtGravitationalParameterOfCentralBody, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        default:
            break;
        }
    }
    // Create partial function if parameter is PPN parameter beta or gamma
    else if( parameter->getParameterName( ).second.first == "global_metric"  )
    {
        switch( parameter->getParameterName( ).first )
        {
        case estimatable_parameters::ppn_parameter_gamma:
            partialFunction = std::bind( &RelativisticAccelerationPartial::wrtPpnParameterGamma, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        case estimatable_parameters::ppn_parameter_beta:
            partialFunction = std::bind( &RelativisticAccelerationPartial::wrtPpnParameterBeta, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        case estimatable_parameters::ppn_parameter_alpha1:
            partialFunction = std::bind( &RelativisticAccelerationPartial::wrtPpnParameterAlpha1, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        case estimatable_parameters::ppn_parameter_alpha2:
            partialFunction = std::bind( &RelativisticAccelerationPartial::wrtPpnParameterAlpha2, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        default:
            break;
        }
    }
    return std::make_pair( partialFunction, numberOfRows );
}

//! Function for updating partial w.r.t. the bodies' states
void RelativisticAccelerationPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        currentRelativeState_ = ( acceleratedBodyState_( ) - centralBodyState_( ) );
        currentSchwarzschildAcceleration_ = currentSchwarzschildAccelerationFunction_( );
        currentSchwarzschildAlphaTermsAcceleration_ = currentSchwarzschildAlphaTermsAccelerationFunction_( );

        computePartialOfSchwarzschildAccelerationCorrectionWrtPosition(
                    currentRelativeState_, currentSchwarzschildAcceleration_, currentSchwarzschildPartialWrtPosition_,
                    centralBodyGravitationalParameterFunction_( ),
                    ppnGammaParameterFunction_( ), ppnBetaParameterFunction_( ) );
        computePartialOfSchwarzschildAccelerationCorrectionWrtVelocity(
                    currentRelativeState_, currentSchwarzschildPartialWrtVelocity_,
                    centralBodyGravitationalParameterFunction_( ),
                    ppnGammaParameterFunction_( ) );


        computePartialOfSchwarzschildAlphaTermsWrtPosition(
                    currentRelativeState_, currentSchwarzschildAlphaTermsAcceleration_,
                    currentSchwarzschildAlphaTermsPartialWrtPosition_,
                    centralBodyGravitationalParameterFunction_( ), acceleratedBodyGravitationalParameterFunction_( ),
                    ppnAlpha1ParameterFunction_( ), ppnAlpha2ParameterFunction_( ));
        computePartialOfSchwarzschildAlphaTermsWrtVelocity(
                    currentRelativeState_, currentSchwarzschildAlphaTermsAcceleration_,
                    currentSchwarzschildAlphaTermsPartialWrtVelocity_,
                    centralBodyGravitationalParameterFunction_( ), acceleratedBodyGravitationalParameterFunction_( ),
                    ppnAlpha1ParameterFunction_( ), ppnAlpha2ParameterFunction_( ));


        if (calculateLenseThirringCorrection_){

            currentLenseThirringAcceleration_ = currentLenseThirringAccelerationFunction_( );
            centralBodyAngularMomentum_ = centralBodyAngularMomentumFunction_( );
//            std::cout<<centralBodyAngularMomentum_<<std::endl;

            computePartialOfLenseThirringAccelerationCorrectionWrtPosition(
                        currentRelativeState_, currentLenseThirringAcceleration_, centralBodyAngularMomentum_,
                        currentLenseThirringPartialWrtPosition_, centralBodyGravitationalParameterFunction_( ),
                        ppnGammaParameterFunction_( ) );
            computePartialOfLenseThirringAccelerationCorrectionWrtVelocity(
                        currentRelativeState_, centralBodyAngularMomentum_,
                        currentLenseThirringPartialWrtVelocity_, centralBodyGravitationalParameterFunction_( ),
                        ppnGammaParameterFunction_( ) );
        } else{
            currentLenseThirringPartialWrtPosition_ = Eigen::Matrix3d::Zero();
            currentLenseThirringPartialWrtVelocity_ = Eigen::Matrix3d::Zero();
        }


//        std::cout<<"pos: "
//                <<currentSchwarzschildPartialWrtPosition_.transpose()
//                <<" // "<<currentSchwarzschildAlphaTermsPartialWrtPosition_.transpose()
//                <<" // "<<currentLenseThirringPartialWrtPosition_.transpose()<<std::endl;

//        std::cout<<"vel: "
//                <<currentSchwarzschildPartialWrtVelocity_.transpose()
//                <<" // "<<currentSchwarzschildAlphaTermsPartialWrtVelocity_.transpose()
//                <<" // "<<currentLenseThirringPartialWrtVelocity_.transpose()<<std::endl;

        currentPartialWrtPosition_ =
                currentSchwarzschildPartialWrtPosition_
                + currentSchwarzschildAlphaTermsPartialWrtPosition_
                + currentLenseThirringPartialWrtPosition_;
        currentPartialWrtVelocity_ =
                currentSchwarzschildPartialWrtVelocity_
                + currentSchwarzschildAlphaTermsPartialWrtVelocity_
                + currentLenseThirringPartialWrtVelocity_;

//        std::cout<<"position partial diagonal: "<<
//                   currentSchwarzschildAlphaTermsPartialWrtPosition_.diagonal().norm()/currentPartialWrtPosition_.diagonal().norm()<<std::endl;
//        std::cout<<"velocity partial diagonal: "<<
//                   currentSchwarzschildAlphaTermsPartialWrtVelocity_.diagonal().norm()/currentPartialWrtVelocity_.diagonal().norm()<<std::endl;

//        std::cout<<"position partial diagonal: "<<
//                   currentLenseThirringPartialWrtPosition_.diagonal().norm()/currentPartialWrtPosition_.diagonal().norm()<<std::endl;
//        std::cout<<"velocity partial diagonal: "<<
//                   currentLenseThirringPartialWrtVelocity_.diagonal().norm()/currentPartialWrtVelocity_.diagonal().norm()<<std::endl;

        currentTime_ = currentTime;
    }
}

}

}


