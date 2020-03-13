#ifndef TIMEVARYINGGRAVITATIONALPARAMETERPARTIAL_H
#define TIMEVARYINGGRAVITATIONALPARAMETERPARTIAL_H

#include "tudatApplications/thesis/MyApplications/timeVaryingGravitationalParameterAcceleration.h"


#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

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
        const double currentTime );


//! Function to compute partial of TVGP w.r.t. gravitational parameter of central body
void computePartialOfTVGPWrtGravitationalParameter(
        const Eigen::Vector6d& currentRelativeState,
        Eigen::MatrixXd& partialMatrix,
        const double gravitationalParameter,
        const double timeVaryingGravitationalParameter,
        const double currentTime );


//! Function to compute partial of TVGP w.r.t. time varying gravitational parameter itself
void computePartialOfTVGPWrtTimeVaryingGravitationalParameter(
        const Eigen::Vector6d& currentRelativeState,
        Eigen::MatrixXd& partialMatrix,
        const double gravitationalParameter,
        const double timeVaryingGravitationalParameter,
        const double currentTime );





//! Class to calculate the partials of the time varying gravitational parameter accelerationn w.r.t. parameters and states.
class TimeVaryingGravitationalParameterPartial: public AccelerationPartial
{
public:

    //! Constructor.
    TimeVaryingGravitationalParameterPartial(
            const std::shared_ptr< TimeVaryingGravitationalParameterAcceleration > accelerationModel,
            const std::string& acceleratedBody,
            const std::string& acceleratingBody):
        AccelerationPartial( acceleratedBody,
                             acceleratingBody,
                             basic_astrodynamics::time_varying_gravitational_parameter_acceleration )
    {
        centralBodyState_  = accelerationModel->getStateFunctionOfCentralBody( );
        acceleratedBodyState_ = accelerationModel->getStateFunctionOfAcceleratedBody( );
        centralBodyGravitationalParameterFunction_ = accelerationModel->getGravitationalParameterFunctionOfCentralBody( );
        timeVaryingGravitationalParameterFunction_ = accelerationModel->getTimeVaryingGravitationalParameterFunction( );

        currentAccelerationFunction_ = std::bind( &TimeVaryingGravitationalParameterAcceleration::getAcceleration,
                                                    accelerationModel );
    }

    //! Destructor
    ~TimeVaryingGravitationalParameterPartial(){};


    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration.
    void wrtPositionOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
    }


    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration.
    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
    }


    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );


    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }


    //! Function to compute partial derivative w.r.t. central body gravitational parameter.
    void wrtGravitationalParameterOfCentralBody( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfTVGPWrtGravitationalParameter(
                    currentRelativeState_,
                    partialMatrix,
                    centralBodyGravitationalParameterFunction_( ),
                    timeVaryingGravitationalParameterFunction_( ),
                    currentTime_);
    }

    //! Function to compute partial derivative w.r.t. time varying gravitational parameter.
    void wrtTimeVaryingGravitationalParameter( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfTVGPWrtGravitationalParameter(
                    currentRelativeState_,
                    partialMatrix,
                    centralBodyGravitationalParameterFunction_( ),
                    timeVaryingGravitationalParameterFunction_( ),
                    currentTime_);
    }

    //! Function for updating partial w.r.t. the bodies' states
    void update( const double currentTime = TUDAT_NAN );


private:

    //! Function to retrieve current state of body exerting acceleration.
    std::function< Eigen::Vector6d( ) > centralBodyState_;

    //! Function to retrieve current state of body undergoing acceleration.
    std::function< Eigen::Vector6d( ) > acceleratedBodyState_;

    //! Function to retrieve current gravitational parameter of central body.
    std::function< double( ) > centralBodyGravitationalParameterFunction_;

    //! Function to retrieve time varying gravitational parameter.
    std::function< double( ) > timeVaryingGravitationalParameterFunction_;

    //! Current partial of relativistic acceleration correction w.r.t. position of body undergoing acceleration
    /*!
     *  Current partial of relativistic acceleration correction w.r.t. position of body undergoing acceleration
     * ( = -partial of central gravity acceleration w.r.t. position of body exerting acceleration),
     *  calculated and set by update( ) function.
     */
    Eigen::Matrix3d currentPartialWrtPosition_;

    //! Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
    Eigen::Vector6d currentRelativeState_;

    //! Current relativistic acceleration correction
    Eigen::Vector3d currentAcceleration_;

    //! Function to retrieve current relativistic acceleration correction.
    std::function< Eigen::Vector3d( ) > currentAccelerationFunction_;



};

} // namespace acceleration_partials

} // namespace tudat

#endif // TIMEVARYINGGRAVITATIONALPARAMETERPARTIAL_H
