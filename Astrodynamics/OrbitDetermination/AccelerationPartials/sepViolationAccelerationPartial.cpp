#include "sepViolationAccelerationPartial.h"


namespace tudat
{
namespace acceleration_partials
{



//! Function to compute partial of TVGP w.r.t. position of body undergoing acceleration
void computePartialOfSEPViolationAccelerationWrtPosition(
        const Eigen::Vector3d& relativePosition,
        const Eigen::Vector3d& sepPositionCorrection,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter)
{


    // if Nordtvedt parameter = 0, there will be no correction, and the corrected position will be equal to the original position
    // as a result the partial matrix will be the zero matrix (terms 1-2 cancel terms 4-5, and term 3 contains a cross product of an identical vector = 0)
    if (sepPositionCorrection.norm( ) != 0.0){
        std::runtime_error( "ERROR: if nordtvedt parameter = 0 exactly, no partial w.r.t. the state can be calculated" );
    }

//    std::cout<<relativePosition<<std::endl;
//    std::cout<<sepPositionCorrection<<std::endl;

    Eigen::Vector3d relativeCorrectedPosition = relativePosition + sepPositionCorrection;
    double relativeNorm = relativePosition.norm( );
    double invSquareRelativeNorm = 1.0 / (relativeNorm * relativeNorm);
    double invCubedRelativeNorm = invSquareRelativeNorm / relativeNorm;

    double relativeCorrectedNorm = relativeCorrectedPosition.norm( );
    double invSquareRelativeCorrectedNorm = 1.0 / (relativeCorrectedNorm * relativeCorrectedNorm);
    double invCubedRelativeCorrectedNorm = invSquareRelativeCorrectedNorm / relativeCorrectedNorm;       


    partialMatrix = Eigen::Matrix3d::Zero( );

    partialMatrix += -Eigen::Matrix3d::Identity( ) * invCubedRelativeCorrectedNorm;
    partialMatrix += 3.0 * invSquareRelativeCorrectedNorm * invCubedRelativeCorrectedNorm
            * relativeCorrectedPosition * relativePosition.transpose( );

    partialMatrix += 3.0 * invSquareRelativeCorrectedNorm * invCubedRelativeCorrectedNorm
            * relativeCorrectedPosition * sepPositionCorrection.transpose( );

    partialMatrix += Eigen::Matrix3d::Identity( ) * invCubedRelativeNorm;
    partialMatrix += -3.0 * invSquareRelativeNorm * invCubedRelativeNorm
            * relativePosition * relativePosition.transpose( );

    partialMatrix *= gravitationalParameter;

//    std::cout<<partialMatrix<<std::endl;


//    partialMatrix = gravitationalParameter *
//            ( Eigen::Matrix3d::Identity( ) * invCubeOfPositionNorm -
//              ( 3.0 * invSquareOfPositionNorm * invCubeOfPositionNorm ) * relativePosition * relativePosition.transpose( ) );

}




//! Function to compute partial of TVGP w.r.t. gravitational parameter of central body
void computePartialOfSEPViolationAccelerationWrtGravitationalParameter(
        const Eigen::Vector3d& relativePosition,
        Eigen::MatrixXd& partialMatrix,
        const Eigen::Vector3d& sepCorrectedRelativePosition)
{
    double positionNorm = relativePosition.norm( );
    double positionNormSEPCorrected = sepCorrectedRelativePosition.norm( );

    //first term, from acceleration due to SEP violation
    partialMatrix = sepCorrectedRelativePosition / (positionNormSEPCorrected
                                                         *positionNormSEPCorrected
                                                         *positionNormSEPCorrected);

    //second term, from conventional acceleration
    partialMatrix += -relativePosition / ( positionNorm * positionNorm * positionNorm );


}




//! Function to compute partial of TVGP w.r.t. the Nordtvedt parameter
void computePartialOfSEPViolationAccelerationWrtNordtvedtParameter(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix)
{
    partialMatrix = nordtvedtPartial;
};

//! Function to compute partial of TVGP w.r.t. the Nordtvedt parameter
void computePartialOfSEPViolationAccelerationWrtPpnParameterGamma(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix)
{
    partialMatrix = -1.0*nordtvedtPartial;
};

//! Function to compute partial of TVGP w.r.t. the Nordtvedt parameter
void computePartialOfSEPViolationAccelerationWrtPpnParameterBeta(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix)
{
    partialMatrix = 4.0*nordtvedtPartial;
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
//        case estimatable_parameters::ppn_parameter_gamma:
//            partialFunction = std::bind( &SEPViolationAccelerationPartial::wrtPpnParameterGamma, this, std::placeholders::_1 );
//            numberOfRows = 1;
//            break;
//        case estimatable_parameters::ppn_parameter_beta:
//            partialFunction = std::bind( &SEPViolationAccelerationPartial::wrtPpnParameterBeta, this, std::placeholders::_1 );
//            numberOfRows = 1;
//            break;
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
