#include "timeVaryingGravitationalParameterPartial.h"



namespace tudat
{

namespace acceleration_partials
{



//! Function to compute partial of TVGP w.r.t. position of body undergoing acceleration
void computePartialOfTVGPWrtPosition(
        const Eigen::Vector6d& currentRelativeState,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter,
        const double timeVaryingGravitationalParameter,
        const double currentTime )
{
    Eigen::Vector3d position = currentRelativeState.segment( 0, 3 );
    double distance = position.norm( );
    double distanceSquared = distance*distance;
    double distanceCubed = distance*distance*distance;

    partialMatrix = gravitationalParameter
            * timeVaryingGravitationalParameter
            * (currentTime/tudat::physical_constants::JULIAN_YEAR)
            * ( Eigen::Matrix3d::Identity( ) / (distanceCubed)
                + 3.0 * (1/distanceSquared) * (1/distanceCubed) * position * position.transpose()
               );

}


//! Function to compute partial of TVGP w.r.t. gravitational parameter of central body
void computePartialOfTVGPWrtGravitationalParameter(
        const Eigen::Vector6d& currentRelativeState,
        Eigen::MatrixXd& partialMatrix,
        const double timeVaryingGravitationalParameter,
        const double currentTime )
{

    Eigen::Vector3d position = currentRelativeState.segment( 0, 3 );
    double distance = position.norm( );
    double distanceCubed = distance*distance*distance;

    partialMatrix = timeVaryingGravitationalParameter
            * (currentTime/tudat::physical_constants::JULIAN_YEAR)
            * position
            / distanceCubed;
}


//! Function to compute partial of TVGP w.r.t. time varying gravitational parameter itself
void computePartialOfTVGPWrtTimeVaryingGravitationalParameter(
        const Eigen::Vector6d& currentRelativeState,
        Eigen::MatrixXd& partialMatrix,
        const double gravitationalParameter,
        const double currentTime )
{

    Eigen::Vector3d position = currentRelativeState.segment( 0, 3 );
    double distance = position.norm( );
    double distanceCubed = distance*distance*distance;

    partialMatrix = gravitationalParameter
            * (currentTime/tudat::physical_constants::JULIAN_YEAR)
            * position
            / distanceCubed;
}




//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
TimeVaryingGravitationalParameterPartial::getParameterPartialFunction(
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
            partialFunction = std::bind( &TimeVaryingGravitationalParameterPartial::wrtGravitationalParameterOfCentralBody, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        case estimatable_parameters::time_varying_gravitational_parameter:
            partialFunction = std::bind( &TimeVaryingGravitationalParameterPartial::wrtTimeVaryingGravitationalParameter, this, std::placeholders::_1);
            numberOfRows = 1;
            break;
        default:
            break;
        }
    }


    return std::make_pair( partialFunction, numberOfRows );
}



// verplaatst naar update() in .h file

// //! Function for updating partial w.r.t. the bodies' states
//void TimeVaryingGravitationalParameterPartial::update( const double currentTime ) // 13/mar: Hier gaat het nu mis, functies toevoegen in constructor?
//{
//    if( !( currentTime_ == currentTime ) )
//    {
//        currentRelativeState_ = ( acceleratedBodyState_( ) - centralBodyState_( ) );
// //        currentAcceleration_ = currentAccelerationFunction_( );

//        currentTime_ = currentTime;

//        computePartialOfTVGPWrtPosition(
//                    currentRelativeState_,
//                    currentPartialWrtPosition_,
//                    centralBodyGravitationalParameterFunction_( ),
//                    timeVaryingGravitationalParameterFunction_( ),
//                    currentTime_);

//    }
//}

}

}
