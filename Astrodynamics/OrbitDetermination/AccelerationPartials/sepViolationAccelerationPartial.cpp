#include "sepViolationAccelerationPartial.h"


namespace tudat
{
namespace acceleration_partials
{



//! Function to compute partial of TVGP w.r.t. position of body undergoing acceleration
void computePartialOfSEPViolationAccelerationWrtPosition(
        const Eigen::Vector3d& relativePosition,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter)
{

    double relativePositionNorm = relativePosition.norm( );
    double invSquareOfPositionNorm = 1.0 / ( relativePositionNorm * relativePositionNorm );
    double invCubeOfPositionNorm = invSquareOfPositionNorm / relativePositionNorm;
    partialMatrix = -gravitationalParameter *
            ( Eigen::Matrix3d::Identity( ) * invCubeOfPositionNorm -
              ( 3.0 * invSquareOfPositionNorm * invCubeOfPositionNorm ) * relativePosition * relativePosition.transpose( ) );

}




//! Function to compute partial of TVGP w.r.t. gravitational parameter of central body
void computePartialOfSEPViolationAccelerationWrtGravitationalParameter(
        const Eigen::Vector3d& relativePosition,
        Eigen::MatrixXd& partialMatrix,
        const Eigen::Vector3d& sepCorrectedRelativePosition)
{
    double positionNorm = relativePosition.norm( );
    double positionNormSEPCorrected = sepCorrectedRelativePosition.norm( );

    //first term, from conventional acceleration
    partialMatrix = relativePosition / ( positionNorm * positionNorm * positionNorm );

    //second term, from acceleration due to SEP violation
    partialMatrix += 0.5*sepCorrectedRelativePosition / (positionNormSEPCorrected*positionNormSEPCorrected);
}




//! Function to compute partial of TVGP w.r.t. the Nordtvedt parameter
void computePartialOfSEPViolationAccelerationWrtNordtvedtParameter(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix)
{
    partialMatrix = nordtvedtPartial;
};




//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
SEPViolationAccelerationPartial::getParameterPartialFunction(
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
            partialFunction = std::bind( &SEPViolationAccelerationPartial::wrtGravitationalParameterOfCentralBody, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        default:
            break;
        }
    }
    // Create partial function if parameter is the Nordtvedt parameter
    else if( parameter->getParameterName( ).second.first == "global_metric"  )
    {
        switch( parameter->getParameterName( ).first )
        {
        case estimatable_parameters::ppn_nordtvedt_parameter:
            partialFunction = std::bind( &SEPViolationAccelerationPartial::wrtNordtvedtParameter, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        default:
            break;
        }
    }
    return std::make_pair( partialFunction, numberOfRows );
}



}
}
