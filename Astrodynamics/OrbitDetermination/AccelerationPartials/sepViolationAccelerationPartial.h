#ifndef SEPVIOLATIONACCELERATIONPARTIAL_H
#define SEPVIOLATIONACCELERATIONPARTIAL_H


#include <memory>
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "tudatApplications/thesis/MyApplications/sepViolationAcceleration.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

namespace tudat
{
namespace acceleration_partials
{



//! Function to compute partial of TVGP w.r.t. position of body undergoing acceleration
void computePartialOfSEPViolationAccelerationWrtPosition(
        const Eigen::Vector3d& relativePosition,
        const Eigen::Vector3d& sepPositionCorrection,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter);



//! Function to compute partial of TVGP w.r.t. gravitational parameter of central body
void computePartialOfSEPViolationAccelerationWrtGravitationalParameter(
        const Eigen::Vector3d& centralBodyPositionShort,
        Eigen::MatrixXd& partialMatrix,
        const Eigen::Vector3d& acceleratedBodyPositionShort,
        const Eigen::Vector3d& sepCorrectionShort,
        const double gravitationalParameterOfCentralBodyShort,
        const Eigen::Vector3d& currentAccelerationShort);


//! Function to compute partial of SEP violation acceleration w.r.t. the Nordtvedt parameter
void computePartialOfSEPViolationAccelerationWrtNordtvedtParameter(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix);

//! Function to compute partial of SEP violation acceleration w.r.t. PPN parameter gamma
void computePartialOfSEPViolationAccelerationWrtPpnParameterGamma(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix);

//! Function to compute partial of SEP violation acceleration w.r.t. PPN parameter beta
void computePartialOfSEPViolationAccelerationWrtPpnParameterBeta(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix);

//! Function to compute partial of SEP violation acceleration w.r.t. PPN parameter alpha1
void computePartialOfSEPViolationAccelerationWrtPpnParameterAlpha1(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix);

//! Function to compute partial of SEP violation acceleration w.r.t. PPN parameter alpha2
void computePartialOfSEPViolationAccelerationWrtPpnParameterAlpha2(
        const Eigen::Vector3d& nordtvedtPartial,
        Eigen::MatrixXd& partialMatrix);


//! Function to compute partial of SEP violation acceleration w.r.t. the Nordtvedt parameter,
//! which is used as a base for the functions above
void calculateNordtvedtPartial(
        const Eigen::Vector3d& correctedPosition,
        Eigen::Vector3d& partialMatrix,
        const Eigen::Vector3d& sepCorrection,
        const double gravitationalParameterOfCentralBody,
        const double nordtvedtParameter);




class SEPViolationAccelerationPartial: public AccelerationPartial
{
public:

    //! Constructor
    SEPViolationAccelerationPartial(
            const std::shared_ptr< relativity::SEPViolationAcceleration > accelerationModel,
            const std::string& acceleratedBody,
            const std::string& acceleratingBody):
        AccelerationPartial( acceleratedBody,
                             acceleratingBody,
                             basic_astrodynamics::sep_violation_acceleration )
    {
        centralBodyPosition_ = accelerationModel->getPositionFunctionOfCentralBody( );
        acceleratedBodyPosition_ = accelerationModel->getPositionFunctionOfAcceleratedBody( );

        sepPositionCorrection_ = accelerationModel->getSEPPositionCorrectionFunction( );

        useNordtvedtConstraint_ = accelerationModel->getUseNordtvedtConstraintFunction( );

        centralBodyGravitationalParameterFunction_ = accelerationModel->getGravitationalParameterFunctionOfCentralBody( );

        currentAcceleration_ = accelerationModel->getAcceleration( );

        nordtvedtParameterFunction_ = accelerationModel->getNordtvedtParameterFunction( );

        nordtvedtPartialFunction_ = accelerationModel->getNordtvedtPartialFunction( );


    }

    //! Destructor
    ~SEPViolationAccelerationPartial(){};


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



    //! Function for determining if the acceleration is dependent on a non-translational integrated state.
    /*!
     *  Function for determining if the acceleration is dependent on a non-translational integrated state.
     *  No dependency is implemented, but a warning is provided if partial w.r.t. mass of body exerting acceleration
     *  (and undergoing acceleration if mutual attraction is used) is requested.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        if( ( ( stateReferencePoint.first == acceleratingBody_ ||
                ( stateReferencePoint.first == acceleratedBody_  ) )
              && integratedStateType == propagators::body_mass_state ) )
        {
            throw std::runtime_error( "Warning, dependency of relativistic acceleration on body masses not yet implemented" );
        }
        return 0;
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
        computePartialOfSEPViolationAccelerationWrtGravitationalParameter(
                    centralBodyPosition_( ),
                    partialMatrix,
                    acceleratedBodyPosition_( ),
                    sepCorrection_,
                    centralBodyGravitationalParameter_,
                    currentAcceleration_);
    }

    //! Function to compute partial derivative w.r.t. nordtvedtparameter.
    void wrtNordtvedtParameter( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfSEPViolationAccelerationWrtNordtvedtParameter(
                    currentNordtvedtPartial_,
                    partialMatrix);
    }

    void wrtPpnParameterGamma( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfSEPViolationAccelerationWrtPpnParameterGamma(
                    currentNordtvedtPartial_,
                    partialMatrix);
    }

    void wrtPpnParameterBeta( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfSEPViolationAccelerationWrtPpnParameterBeta(
                    currentNordtvedtPartial_,
                    partialMatrix);
    }

    void wrtPpnParameterAlpha1( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfSEPViolationAccelerationWrtPpnParameterAlpha1(
                    currentNordtvedtPartial_,
                    partialMatrix);
    }

    void wrtPpnParameterAlpha2( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfSEPViolationAccelerationWrtPpnParameterAlpha2(
                    currentNordtvedtPartial_,
                    partialMatrix);
    }



    //! Function for updating partial w.r.t. the bodies' states
    void update( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime ) )
        {
            sepCorrection_ = sepPositionCorrection_( ) ;

            sepCorrectedCentralBodyPosition_ =
                    centralBodyPosition_( ) + sepCorrection_ ;

            centralBodyGravitationalParameter_ =
                    centralBodyGravitationalParameterFunction_( );

            nordtvedtParameter_ = nordtvedtParameterFunction_( );

            currentRelativePosition_ =
                    ( acceleratedBodyPosition_( ) - centralBodyPosition_( ) );
            currentSEPCorrectedRelativePosition_ =
                    ( acceleratedBodyPosition_( ) - sepCorrectedCentralBodyPosition_ );



            computePartialOfSEPViolationAccelerationWrtPosition(
                        currentRelativePosition_,
                        sepCorrection_,
                        currentPartialWrtPosition_,
                        centralBodyGravitationalParameter_
                        );

//            currentNordtvedtPartial_ = nordtvedtPartialFunction_( ); // uses the partial as given in Genova et al 2018, Nature communications, eq. 10

            calculateNordtvedtPartial( // uses a self-derived partial equation based on the simplified delta r_SEP
                        currentSEPCorrectedRelativePosition_,
                        currentNordtvedtPartial_,
                        sepCorrection_,
                        centralBodyGravitationalParameter_,
                        nordtvedtParameter_);

        }
    }





private:

    //! Function to retrieve current state of body exerting acceleration.
    std::function< Eigen::Vector3d( ) > centralBodyPosition_;

    //! Function to retrieve current state of body undergoing acceleration.
    std::function< Eigen::Vector3d( ) > acceleratedBodyPosition_;

    //! Function to retrieve current state of body exerting acceleration.
    std::function< Eigen::Vector3d( ) > sepPositionCorrection_;

    //! Function to retrieve current state of body undergoing acceleration.
    std::function< Eigen::Vector3d( ) > nordtvedtPartialFunction_;

    std::function< bool( ) > useNordtvedtConstraint_;


    //! Function to retrieve current gravitational parameter of central body.
    std::function< double( ) > centralBodyGravitationalParameterFunction_;

    std::function< double( ) > nordtvedtParameterFunction_;

    //! Current partial w.r.t. position of body undergoing acceleration
    Eigen::Matrix3d currentPartialWrtPosition_;

    //! Relative position of body undergoing acceleration.
    Eigen::Vector3d currentRelativePosition_;

    //! Relative position of body undergoing acceleration.
    Eigen::Vector3d currentSEPCorrectedRelativePosition_;

    double centralBodyGravitationalParameter_;

    double nordtvedtParameter_;

    Eigen::Vector3d sepCorrection_;

    Eigen::Vector3d currentAcceleration_;

    //! Function to retrieve current state of body exerting acceleration.
    Eigen::Vector3d sepCorrectedCentralBodyPosition_;

    //! Relative position of body undergoing acceleration.
    Eigen::Vector3d currentNordtvedtPartial_;




};


}
}

#endif // SEPVIOLATIONACCELERATIONPARTIAL_H
