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
    // as a result the partial matrix will be the zero matrix (the terms with and without correction are then equal, but with opposite signs)
    if (sepPositionCorrection.norm( ) == 0.0){
        partialMatrix = Eigen::Matrix3d::Zero();
    } else{

        // convert doubles to long doubles
        Eigen::Matrix<long double, 3, 1> relativePositionLong = relativePosition.cast<long double>();
        Eigen::Matrix<long double, 3, 1> sepPositionCorrectionLong = sepPositionCorrection.cast<long double>();
        long double gravitationalParameterLong = static_cast<long double>(gravitationalParameter);

        // get necessary variables
        Eigen::Matrix<long double, 3, 1> relativeCorrectedPosition = relativePositionLong - sepPositionCorrectionLong;
        long double relativeNorm = relativePositionLong.norm( );
        long double invSquareRelativeNorm = 1.0 / (relativeNorm * relativeNorm);
        long double invCubedRelativeNorm = invSquareRelativeNorm / relativeNorm;

        long double relativeCorrectedNorm = relativeCorrectedPosition.norm( );
        long double invSquareRelativeCorrectedNorm = 1.0 / (relativeCorrectedNorm * relativeCorrectedNorm);
        long double invCubedRelativeCorrectedNorm = invSquareRelativeCorrectedNorm / relativeCorrectedNorm;

        // calculate partial derivative terms
        Eigen::Matrix<long double, 3, 3> partialMatrixLong = Eigen::Matrix<long double, 3, 3>::Zero( );

        partialMatrixLong += -1.0 * Eigen::Matrix<long double, 3, 3>::Identity( ) * invCubedRelativeCorrectedNorm;
        partialMatrixLong += 3.0 * invSquareRelativeCorrectedNorm * invCubedRelativeCorrectedNorm
                * relativeCorrectedPosition * relativeCorrectedPosition.transpose( );

        partialMatrixLong += Eigen::Matrix<long double, 3, 3>::Identity( ) * invCubedRelativeNorm;
        partialMatrixLong += -3.0 * invSquareRelativeNorm * invCubedRelativeNorm
                * relativePositionLong * relativePositionLong.transpose( );

        partialMatrixLong *= gravitationalParameterLong;

        // cast back to double for output
        partialMatrix = partialMatrixLong.cast<double>();

    }
}




//! Function to compute partial of TVGP w.r.t. gravitational parameter of central body
void computePartialOfSEPViolationAccelerationWrtGravitationalParameter(
        const Eigen::Vector3d& centralBodyPositionShort,
        Eigen::MatrixXd& partialMatrix,
        const Eigen::Vector3d& acceleratedBodyPositionShort,
        const Eigen::Vector3d& sepCorrectionShort,
        const double gravitationalParameterOfCentralBodyShort,
        const Eigen::Vector3d& currentAccelerationShort)
{

    // convert everything to long for additional accuracy
    Eigen::Matrix<long double, 3, 1> centralBodyPosition = centralBodyPositionShort.cast<long double>();
    Eigen::Matrix<long double, 3, 1> acceleratedBodyPosition = acceleratedBodyPositionShort.cast<long double>();
    Eigen::Matrix<long double, 3, 1> sepCorrection = sepCorrectionShort.cast<long double>();
    long double gravitationalParameterOfCentralBody = static_cast<long double>(gravitationalParameterOfCentralBodyShort);
    Eigen::Matrix<long double, 3, 1> currentAcceleration = currentAccelerationShort.cast<long double>();

//    partialMatrix = Eigen::Vector3d::Zero();

    // some useful variables
    Eigen::Matrix<long double, 3, 1> relativePosition =
            acceleratedBodyPosition - centralBodyPosition;
    long double positionNorm = relativePosition.norm();
    long double invCubedRelativeNorm = 1.0 / (positionNorm * positionNorm * positionNorm);

    Eigen::Matrix<long double, 3, 1> correctedCentralBodyPosition =
            centralBodyPosition + sepCorrection;
    Eigen::Matrix<long double, 3, 1> sepCorrectedRelativePosition =
            acceleratedBodyPosition - correctedCentralBodyPosition;
    long double positionNormSEPCorrected = sepCorrectedRelativePosition.norm( );
    long double invSquareRelativeCorrectedNorm = 1.0 / (positionNormSEPCorrected * positionNormSEPCorrected);
    long double invCubedRelativeCorrectedNorm = invSquareRelativeCorrectedNorm / positionNormSEPCorrected;

//    Eigen::Matrix<long double, 3, 1> c = sepCorrection / gravitationalParameterOfCentralBody;

//    Eigen::Matrix<long double, 3, 1> currentConventionalAcceleration
//            = gravitation::computeGravitationalAccelerationLong(
//                acceleratedBodyPosition,
//                gravitationalParameterOfCentralBody,
//                centralBodyPosition);
//    Eigen::Matrix<long double, 3, 1> currentCorrectedAcceleration
//            = gravitation::computeGravitationalAccelerationLong(
//                acceleratedBodyPosition,
//                gravitationalParameterOfCentralBody,
//                correctedCentralBodyPosition);

//    Eigen::Matrix<long double, 3, 1> currentAcceleration =
//            currentCorrectedAcceleration - currentConventionalAcceleration;


    // calculate terms
    Eigen::Matrix<long double, 3, 1> term1 =
            currentAcceleration / gravitationalParameterOfCentralBody;

    Eigen::Matrix<long double, 3, 3> term2a =
            Eigen::Matrix<long double, 3, 3>::Identity( ) * invCubedRelativeCorrectedNorm;

    Eigen::Matrix<long double, 3, 3> term2b =
            3.0 * invSquareRelativeCorrectedNorm * invCubedRelativeCorrectedNorm
            * sepCorrectedRelativePosition * sepCorrectedRelativePosition.transpose( );

    Eigen::Matrix<long double, 3, 1> term2 =
            (term2a - term2b) * sepCorrection;

    Eigen::Matrix<long double, 3, 1> partialMatrixLong = term1 - term2;

    partialMatrix = partialMatrixLong.cast<double>();

//    std::cout<<"partial wrt grav parameter"<<std::endl;
//    std::cout<<gravitationalParameterOfCentralBody<<std::endl;
////    std::cout<<currentConventionalAcceleration<<std::endl;
////    std::cout<<currentCorrectedAcceleration<<std::endl;
//    std::cout<<currentAcceleration<<std::endl;
//    std::cout<<term1<<std::endl;
//    std::cout<<term2a<<std::endl;
//    std::cout<<term2b<<std::endl;
//    std::cout<<term2<<std::endl;
//    std::cout<<partialMatrixLong<<std::endl;
//    std::cout<<partialMatrix<<std::endl;



//    Eigen::Matrix<long double, 3, 1> partialMatrixLong =
//            -1.0 * sepCorrectedRelativePosition * invCubedRelativeCorrectedNorm;

//    partialMatrixLong += 3.0
//            * sepCorrectedRelativePosition.transpose().dot(sepCorrectedRelativePosition)
//            * c * invCubedRelativeCorrectedNorm * invSquareRelativeCorrectedNorm
//            / gravitationalParameterOfCentralBody;

//    partialMatrixLong -= c
//            * invCubedRelativeCorrectedNorm
//            / gravitationalParameterOfCentralBody;

//    partialMatrixLong += relativePosition * invCubedRelativeNorm;


}




//! Function to compute partial of SEP violation acceleration w.r.t. the Nordtvedt parameter
void computePartialOfSEPViolationAccelerationWrtNordtvedtParameter(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix)
{
    partialMatrix = nordtvedtPartial;
//    std::cout<<partialMatrix.transpose()<<std::endl;
};

//! Function to compute partial of SEP violation acceleration w.r.t. PPN parameter gamma
void computePartialOfSEPViolationAccelerationWrtPpnParameterGamma(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix)
{
    partialMatrix = -1.0*nordtvedtPartial;
};

//! Function to compute partial of SEP violation acceleration w.r.t. PPN parameter beta
void computePartialOfSEPViolationAccelerationWrtPpnParameterBeta(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix)
{
    partialMatrix = 4.0*nordtvedtPartial;
};

//! Function to compute partial of SEP violation acceleration w.r.t. PPN parameter alpha1
void computePartialOfSEPViolationAccelerationWrtPpnParameterAlpha1(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix)
{
    partialMatrix = -1.0*nordtvedtPartial;
};

//! Function to compute partial of SEP violation acceleration w.r.t. PPN parameter alpha2
void computePartialOfSEPViolationAccelerationWrtPpnParameterAlpha2(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix)
{
    partialMatrix = (-2.0/3.0)*nordtvedtPartial;
};


//! Function to compute partial of SEP violation acceleration w.r.t. the Nordtvedt parameter,
//! which is used as a base for the functions above
void calculateNordtvedtPartial(
        const Eigen::Vector3d& correctedPosition,
        Eigen::Vector3d& partialMatrix,
        const Eigen::Vector3d& sepCorrection,
        const double gravitationalParameterOfCentralBody,
        const double nordtvedtParameter){

    if (sepCorrection.norm() == 0.0){
        partialMatrix = Eigen::Vector3d::Zero();
    } else{

        double distance = correctedPosition.norm();

        Eigen::Matrix3d bracketTerm = Eigen::Matrix3d::Identity( )
                / (distance * distance * distance);

        bracketTerm -= 3.0 * correctedPosition
                * correctedPosition.transpose()
                / (distance * distance * distance * distance *distance);

        partialMatrix = bracketTerm
                * gravitationalParameterOfCentralBody
                * sepCorrection
                / nordtvedtParameter;
    }
}




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
//        case estimatable_parameters::gravitational_parameter:
//            partialFunction = std::bind( &SEPViolationAccelerationPartial::wrtGravitationalParameterOfCentralBody, this, std::placeholders::_1 );
//            numberOfRows = 1;
//            break;
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
            if (useNordtvedtConstraint_( ) == true){
                partialFunction = std::bind( &SEPViolationAccelerationPartial::wrtPpnParameterGamma, this, std::placeholders::_1 );
                numberOfRows = 1;
            }
            break;
        case estimatable_parameters::ppn_parameter_beta:
            if (useNordtvedtConstraint_( ) == true){
                partialFunction = std::bind( &SEPViolationAccelerationPartial::wrtPpnParameterBeta, this, std::placeholders::_1 );
                numberOfRows = 1;
            }
            break;
        case estimatable_parameters::ppn_parameter_alpha1:
            if (useNordtvedtConstraint_( ) == true){
                partialFunction = std::bind( &SEPViolationAccelerationPartial::wrtPpnParameterAlpha1, this, std::placeholders::_1 );
                numberOfRows = 1;
            }
            break;
        case estimatable_parameters::ppn_parameter_alpha2:
            if (useNordtvedtConstraint_( ) == true){
                partialFunction = std::bind( &SEPViolationAccelerationPartial::wrtPpnParameterAlpha2, this, std::placeholders::_1 );
                numberOfRows = 1;
            }
            break;
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
