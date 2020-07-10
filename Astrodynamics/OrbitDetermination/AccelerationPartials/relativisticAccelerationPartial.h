/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RELATIVISTICACCELERATIONPARTIAL_H
#define TUDAT_RELATIVISTICACCELERATIONPARTIAL_H

#include <memory>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Relativity/relativisticAccelerationCorrection.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function to compute partial of Schwarzschild acceleration correction w.r.t. position of body undergoing acceleration
/*!
 * Function to compute partial of Schwarzschild acceleration correction w.r.t. position of body undergoing acceleration
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param currentAcceleration Current Schwarzschild acceleration correction
 * \param partialMatrix Requested (returnd by reference)
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param ppnParameterGamma PPN parameter gamma
 * \param ppnParameterBeta PPN parameter beta
 */
void computePartialOfSchwarzschildAccelerationCorrectionWrtPosition(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter,
        const double ppnParameterGamma = 1.0,
        const double ppnParameterBeta = 1.0);

//! Function to compute partial of Schwarzschild acceleration correction w.r.t. velocity of body undergoing acceleration
/*!
 * Function to compute partial of Schwarzschild acceleration correction w.r.t. velocity of body undergoing acceleration
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param partialMatrix Requested (returnd by reference)
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param ppnParameterGamma PPN parameter gamma
 */
void computePartialOfSchwarzschildAccelerationCorrectionWrtVelocity(
        const Eigen::Vector6d& relativeState,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter,
        const double ppnParameterGamma = 1.0);

//! Function to compute partial derivative of Schwarzschild acceleration correction w.r.t. central body gravitational patameter.
/*!
 * Function to compute partial derivative of Schwarzschild acceleration correction w.r.t. central body gravitational patameter.
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param partialMatrix Requested (returnd by reference)
 * \param ppnParameterGamma PPN parameter gamma
 * \param ppnParameterBeta PPN parameter bet
 */
void computePartialOfSchwarzschildAccelerationCorrectionWrtGravitationalParameter(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partialMatrix,
        const double ppnParameterGamma = 1.0,
        const double ppnParameterBeta = 1.0 );

//! Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter gamma
/*!
 * Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter gamma
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param partialMatrix Requested (returnd by reference)
 */
void computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterGamma(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partialMatrix );

//! Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter beta
/*!
 * Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter beta
 * \param relativeState Position of body undergoing, w.r.t. body exerting, acceleration.
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param partialMatrix Requested (returnd by reference)
 */
void computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterBeta(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partialMatrix );






//! Function to compute partial of Schwarzschild acceleration correction w.r.t. position of body undergoing acceleration
/*!
 * Function to compute partial of Schwarzschild acceleration correction w.r.t. position of body undergoing acceleration
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param currentAcceleration Current Schwarzschild acceleration correction
 * \param partialMatrix Requested (returnd by reference)
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param ppnParameterGamma PPN parameter gamma
 * \param ppnParameterBeta PPN parameter beta
 */
void computePartialOfSchwarzschildAlphaTermsWrtPosition(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameterBodyExertingAcceleration,
        const double gravitationalParameterBodyUndergoingAcceleration,
        const double ppnParameterAlpha1 = 0.0,
        const double ppnParameterAlpha2 = 0.0);

//! Function to compute partial of Schwarzschild acceleration correction w.r.t. velocity of body undergoing acceleration
/*!
 * Function to compute partial of Schwarzschild acceleration correction w.r.t. velocity of body undergoing acceleration
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param partialMatrix Requested (returnd by reference)
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param ppnParameterGamma PPN parameter gamma
 */
void computePartialOfSchwarzschildAlphaTermsWrtVelocity(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameterBodyExertingAcceleration,
        const double gravitationalParameterBodyUndergoingAcceleration,
        const double ppnParameterAlpha1 = 0.0,
        const double ppnParameterAlpha2 = 0.0);


//! Function to compute partial derivative of Schwarzschild acceleration correction w.r.t. central body gravitational patameter.
/*!
 * Function to compute partial derivative of Schwarzschild acceleration correction w.r.t. central body gravitational patameter.
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param partialMatrix Requested (returnd by reference)
 * \param ppnParameterGamma PPN parameter gamma
 * \param ppnParameterBeta PPN parameter bet
 */
void computePartialOfSchwarzschildAlphaTermsWrtGravitationalParameter(
        Eigen::Vector3d& currentAcceleration,
        Eigen::MatrixXd& partialMatrix,
        const double gravitationalParameterBodyExertingAcceleration);


//! Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter gamma
/*!
 * Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter gamma
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param partialMatrix Requested (returnd by reference)
 */
void computePartialOfSchwarzschildAlphaTermsWrtPpnParameterAlpha1(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        Eigen::MatrixXd& partialMatrix,
        const double gravitationalParameterBodyExertingAcceleration,
        const double gravitationalParameterBodyUndergoingAcceleration);

//! Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter beta
/*!
 * Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter beta
 * \param relativeState Position of body undergoing, w.r.t. body exerting, acceleration.
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param partialMatrix Requested (returnd by reference)
 */
void computePartialOfSchwarzschildAlphaTermsWrtPpnParameterAlpha2(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        Eigen::MatrixXd& partialMatrix,
        const double gravitationalParameterBodyUndergoingAcceleration);








//! Function to compute partial of LenseThirring acceleration correction w.r.t. position of body undergoing acceleration
/*!
 * Function to compute partial of LenseThirring acceleration correction w.r.t. position of body undergoing acceleration
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param currentAcceleration Current Schwarzschild acceleration correction
 * \param partialMatrix Requested (returnd by reference)
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param ppnParameterGamma PPN parameter gamma
 * \param ppnParameterBeta PPN parameter beta
 */
void computePartialOfLenseThirringAccelerationCorrectionWrtPosition(
        const Eigen::Vector6d& relativeState,
        Eigen::Vector3d& currentAcceleration,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter,
        const double ppnParameterGamma = 1.0);

//! Function to compute partial of LenseThirring acceleration correction w.r.t. velocity of body undergoing acceleration
/*!
 * Function to compute partial of LenseThirring acceleration correction w.r.t. velocity of body undergoing acceleration
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param partialMatrix Requested (returnd by reference)
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param ppnParameterGamma PPN parameter gamma
 */
void computePartialOfLenseThirringAccelerationCorrectionWrtVelocity(
        const Eigen::Vector6d& relativeState,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter,
        const double ppnParameterGamma = 1.0 );

//! Function to compute partial derivative of LenseThirring acceleration correction w.r.t. central body gravitational patameter.
/*!
 * Function to compute partial derivative of LenseThirring acceleration correction w.r.t. central body gravitational patameter.
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param gravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param partialMatrix Requested (returnd by reference)
 * \param ppnParameterGamma PPN parameter gamma
 * \param ppnParameterBeta PPN parameter bet
 */
void computePartialOfLenseThirringAccelerationCorrectionWrtGravitationalParameter(
        const Eigen::Vector6d& relativeState,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        Eigen::MatrixXd& partialMatrix,
        const double ppnParameterGamma = 1.0);


//! Function to compute the partial derivative of LenseThirring acceleration correction w.r.t. angular momentum of central body
void computePartialOfLenseThirringAccelerationCorrectionWrtAngularMomentum(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        Eigen::MatrixXd& partialMatrix,
        const double ppnParameterGamma = 1.0);



//! Function to compute the partial derivative of LenseThirring acceleration correction w.r.t. PPN parameter gamma
void computePartialOfLenseThirringAccelerationCorrectionWrtPpnParameterGamma(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        Eigen::MatrixXd& partialMatrix );


// Beta does not appear in the lense-thirring correction



//! Class to calculate the partials of the relativistic acceleration correction w.r.t. parameters and states.
class RelativisticAccelerationPartial: public AccelerationPartial
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param accelerationModel Relativistic acceleration correction w.r.t. which partials are to be taken.
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     */
    RelativisticAccelerationPartial(
            const std::shared_ptr< relativity::RelativisticAccelerationCorrection > accelerationModel,
            const std::string& acceleratedBody,
            const std::string& acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::relativistic_correction_acceleration )
    {
        // Only Schwarzschild correction is implemented so far.
        if( accelerationModel->getCalculateDeSitterCorrection( ) )
        {
            throw std::runtime_error( "Error when creating relativistic acceleration correction partial, de Sitter not yet implemented" );
        }
        centralBodyState_  = accelerationModel->getStateFunctionOfCentralBody( );
        acceleratedBodyState_ = accelerationModel->getStateFunctionOfAcceleratedBody( );

        ppnGammaParameterFunction_ = accelerationModel->getPpnParameterGammaFunction_( );
        ppnBetaParameterFunction_ = accelerationModel->getPpnParameterBetaFunction_( );
        ppnAlpha1ParameterFunction_ = accelerationModel->getPpnParameterAlpha1Function_( );
        ppnAlpha2ParameterFunction_ = accelerationModel->getPpnParameterAlpha2Function_( );

        centralBodyGravitationalParameterFunction_ = accelerationModel->getGravitationalParameterFunctionOfCentralBody( );
        acceleratedBodyGravitationalParameterFunction_ = accelerationModel->getGravitationalParameterFunctionOfAcceleratedBody( );

        currentSchwarzschildAccelerationFunction_ =
                std::bind( &relativity::RelativisticAccelerationCorrection::getSchwarzschildAcceleration,
                           accelerationModel );
        currentSchwarzschildAlphaTermsAccelerationFunction_ =
                std::bind( &relativity::RelativisticAccelerationCorrection::getSchwarzschildAlphaTermsAcceleration,
                           accelerationModel );
        currentLenseThirringAccelerationFunction_ =
                std::bind( &relativity::RelativisticAccelerationCorrection::getLenseThirringAcceleration,
                           accelerationModel );
        centralBodyAngularMomentumFunction_ =
                std::bind( &relativity::RelativisticAccelerationCorrection::getCentralBodyAngularMomentumVector,
                           accelerationModel );

        calculateLenseThirringCorrection_ = accelerationModel->getCalculateLenseThirringCorrection( );
//        centralBodyAngularMomentumFunction_ = accelerationModel->getCentralBodyAngularMomentum_( );
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
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

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
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

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratingBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
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
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Function to compute partial derivative of relativistic acceleration correction w.r.t. central body gravitational patameter.
    /*!
     * Function to compute partial derivative of relativistic acceleration correction w.r.t. central body gravitational patameter.
     * \param partialMatrix Requested (returnd by reference)
     */
    void wrtGravitationalParameterOfCentralBody( Eigen::MatrixXd& partialMatrix )
    {
        Eigen::MatrixXd partialMatrixSchwarzschild;
        Eigen::MatrixXd partialMatrixSchwarzschildAlphaTerms;
        Eigen::MatrixXd partialMatrixLenseThirring;

        computePartialOfSchwarzschildAccelerationCorrectionWrtGravitationalParameter(
                    currentRelativeState_, centralBodyGravitationalParameterFunction_( ),
                    partialMatrixSchwarzschild,
                    ppnGammaParameterFunction_( ), ppnBetaParameterFunction_( ) );

        computePartialOfSchwarzschildAlphaTermsWrtGravitationalParameter(
                currentSchwarzschildAlphaTermsAcceleration_,
                partialMatrixSchwarzschildAlphaTerms,
                centralBodyGravitationalParameterFunction_( ));

        if (calculateLenseThirringCorrection_){
            computePartialOfLenseThirringAccelerationCorrectionWrtGravitationalParameter(
                        currentRelativeState_, centralBodyAngularMomentumFunction_( ),
                        partialMatrixLenseThirring, ppnGammaParameterFunction_( ) );
        } else{
            partialMatrixLenseThirring = Eigen::Vector3d::Zero();
        }

        partialMatrix = partialMatrixSchwarzschild
                + partialMatrixSchwarzschildAlphaTerms
                + partialMatrixLenseThirring;

//        std::cout<<"norm velocity partial: "<<
//                   partialMatrixSchwarzschildAlphaTerms.norm()/partialMatrix.norm()<<std::endl;

//        std::cout<<"wrtGravitationalParameter: "<<partialMatrixSchwarzschild.transpose()
//                <<" // "<<partialMatrixSchwarzschildAlphaTerms.transpose()
//                <<" // "<<partialMatrixLenseThirring.transpose()<<std::endl;
    }


    void wrtAngularMomentumOfCentralBody( Eigen::MatrixXd& partialMatrix )
    {

        computePartialOfLenseThirringAccelerationCorrectionWrtAngularMomentum(
                    currentRelativeState_, centralBodyGravitationalParameterFunction_( ),
                    centralBodyAngularMomentumFunction_( ),
                    partialMatrix, ppnGammaParameterFunction_( ) );

    }


    //! Function to compute partial derivative of relativistic acceleration correction w.r.t. PPN parameter gamma
    /*!
     * Function to compute partial derivative of relativistic acceleration correction w.r.t. PPN parameter gamma
     * \param partialMatrix Requested (returnd by reference)
     */
    void wrtPpnParameterGamma( Eigen::MatrixXd& partialMatrix )
    {
        Eigen::MatrixXd partialMatrixSchwarzschild;
        Eigen::MatrixXd partialMatrixLenseThirring;

        computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterGamma(
               currentRelativeState_,
               centralBodyGravitationalParameterFunction_( ),
               partialMatrixSchwarzschild );

        if (calculateLenseThirringCorrection_){
            computePartialOfLenseThirringAccelerationCorrectionWrtPpnParameterGamma(
                   currentRelativeState_,
                   centralBodyGravitationalParameterFunction_( ),
                   centralBodyAngularMomentumFunction_( ),
                   partialMatrixLenseThirring );
        } else{
            partialMatrixLenseThirring = Eigen::Vector3d::Zero();
        }

        partialMatrix = partialMatrixSchwarzschild
                + partialMatrixLenseThirring;
//        std::cout<<"wrtParameterGamma:"<<partialMatrix.transpose()<<std::endl;
    }

    //! Function to compute partial derivative of relativistic acceleration correction w.r.t. PPN parameter beta
    /*!
     * Function to compute partial derivative of relativistic acceleration correction w.r.t. PPN parameter beta
     * \param partialMatrix Requested (returnd by reference)
     */
    void wrtPpnParameterBeta( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterBeta(
                    currentRelativeState_, centralBodyGravitationalParameterFunction_( ), partialMatrix );
//        std::cout<<"wrtParameterBeta:"<<partialMatrix.transpose()<<std::endl;
        // Parameter beta does not occur in the LenseThirring acceleration correction
    }


    //! Function to compute partial derivative of relativistic acceleration correction w.r.t. PPN parameter gamma
    /*!
     * Function to compute partial derivative of relativistic acceleration correction w.r.t. PPN parameter gamma
     * \param partialMatrix Requested (returnd by reference)
     */
    void wrtPpnParameterAlpha1( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfSchwarzschildAlphaTermsWrtPpnParameterAlpha1(
               currentRelativeState_,
               currentSchwarzschildAlphaTermsAcceleration_,
               partialMatrix,
               centralBodyGravitationalParameterFunction_( ),
               acceleratedBodyGravitationalParameterFunction_( ) );

//        std::cout<<"wrtParameterAlpha1:"<<partialMatrix.transpose()<<std::endl;

    }


    //! Function to compute partial derivative of relativistic acceleration correction w.r.t. PPN parameter gamma
    /*!
     * Function to compute partial derivative of relativistic acceleration correction w.r.t. PPN parameter gamma
     * \param partialMatrix Requested (returnd by reference)
     */
    void wrtPpnParameterAlpha2( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfSchwarzschildAlphaTermsWrtPpnParameterAlpha2(
               currentRelativeState_,
               currentSchwarzschildAlphaTermsAcceleration_,
               partialMatrix,
               acceleratedBodyGravitationalParameterFunction_( ) );

//        std::cout<<"wrtParameterAlpha2:"<<partialMatrix.transpose()<<std::endl;

    }


    Eigen::MatrixXd returnPositionPartialOfOneTerm( const int term)
    {
        Eigen::MatrixXd partial;
        switch(term){
        case 1: { partial = currentSchwarzschildPartialWrtPosition_; break; }
        case 2: { partial = currentSchwarzschildAlphaTermsPartialWrtPosition_; break; }
        case 3: { partial = currentLenseThirringPartialWrtPosition_; break; }
        default: { std::runtime_error("ERROR: invalid input for term"); break; }
        }
        return partial;
    }

    Eigen::MatrixXd returnVelocityPartialOfOneTerm( const int term)
    {
        Eigen::MatrixXd partial;
        switch(term){
        case 1: { partial = currentSchwarzschildPartialWrtVelocity_; break; }
        case 2: { partial = currentSchwarzschildAlphaTermsPartialWrtVelocity_; break; }
        case 3: { partial = currentLenseThirringPartialWrtVelocity_; break; }
        default: { std::runtime_error("ERROR: invalid input for term"); break; }
        }
        return partial;
    }



    //! Function for updating partial w.r.t. the bodies' states
    /*!
     *  Function for updating common blocks of partial to current state. For the this acceleration,
     *  position and velocity partials is computed and set. Also, member variables are updated to current time/state.
     *  \param currentTime Time at which partials are to be calculated
     */
    void update( const double currentTime = TUDAT_NAN );

private:

    //! Function to retrieve current state of body exerting acceleration.
    std::function< Eigen::Vector6d( ) > centralBodyState_;

    //! Function to retrieve current state of body undergoing acceleration.
    std::function< Eigen::Vector6d( ) > acceleratedBodyState_;

    //! Function to retrieve PPN parameter gamma
    std::function< double( ) > ppnGammaParameterFunction_;

    //! Function to retrieve PPN parameter beta
    std::function< double( ) > ppnBetaParameterFunction_;

    //! Function to retrieve PPN parameter alpha1
    std::function< double( ) > ppnAlpha1ParameterFunction_;

    //! Function to retrieve PPN parameter alpha2
    std::function< double( ) > ppnAlpha2ParameterFunction_;

    //! Function to retrieve the angular momentum of the central body
    std::function< Eigen::Vector3d( ) > centralBodyAngularMomentumFunction_;

    //! Function to retrieve current gravitational parameter of central body.
    std::function< double( ) > centralBodyGravitationalParameterFunction_;

    //! Function to retrieve current gravitational parameter of central body.
    std::function< double( ) > acceleratedBodyGravitationalParameterFunction_;

    //! Function to retrieve current relativistic acceleration correction.
    std::function< Eigen::Vector3d( ) > currentSchwarzschildAccelerationFunction_;
    std::function< Eigen::Vector3d( ) > currentSchwarzschildAlphaTermsAccelerationFunction_;
    std::function< Eigen::Vector3d( ) > currentLenseThirringAccelerationFunction_;

    //! Current partial of relativistic acceleration correction w.r.t. position of body undergoing acceleration
    /*!
     *  Current partial of relativistic acceleration correction w.r.t. position of body undergoing acceleration
     * ( = -partial of central gravity acceleration w.r.t. position of body exerting acceleration),
     *  calculated and set by update( ) function.
     */
    Eigen::Matrix3d currentSchwarzschildPartialWrtPosition_;
    Eigen::Matrix3d currentSchwarzschildAlphaTermsPartialWrtPosition_;
    Eigen::Matrix3d currentLenseThirringPartialWrtPosition_;
    Eigen::Matrix3d currentPartialWrtPosition_;

    //! Current partial of relativistic acceleration correction w.r.t. velocity of body undergoing acceleration
    /*!
     *  Current partial of relativistic acceleration correction w.r.t. velocity of body undergoing acceleration
     * ( = -partial of central gravity acceleration w.r.t. velocity of body exerting acceleration),
     *  calculated and set by update( ) function.
     */
    Eigen::Matrix3d currentSchwarzschildPartialWrtVelocity_;
    Eigen::Matrix3d currentSchwarzschildAlphaTermsPartialWrtVelocity_;
    Eigen::Matrix3d currentLenseThirringPartialWrtVelocity_;
    Eigen::Matrix3d currentPartialWrtVelocity_;

    //! Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
    Eigen::Vector6d currentRelativeState_;

    //! Current relativistic acceleration correction
    Eigen::Vector3d currentSchwarzschildAcceleration_;
    Eigen::Vector3d currentSchwarzschildAlphaTermsAcceleration_;
    Eigen::Vector3d currentLenseThirringAcceleration_;

    Eigen::Vector3d centralBodyAngularMomentum_;


    //! Boolean which indicates if Lense Thirring correction is to be included
    bool calculateLenseThirringCorrection_;

};

}

}

#endif // TUDAT_RELATIVISTICACCELERATIONPARTIAL_H
