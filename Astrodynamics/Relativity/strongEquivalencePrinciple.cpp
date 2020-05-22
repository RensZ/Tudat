

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

namespace tudat
{

namespace relativity
{





//! get the gravitational self energy of a body, which is needed in the function below
//! values and equation come from Genova et al. 2018, nature communications
//! radii from //https://nssdc.gsfc.nasa.gov/planetary/factsheet/<INSERTPLANET>fact.html
double getGravitationalSelfEnergy(
        const std::string bodyName,
        const double gravitationalParameter)
{
    double gravitationalSelfEnergy;

    // bodies for which gravitational self energy value is quite well known
    if (bodyName == "Sun"){
        gravitationalSelfEnergy = -3.52E-6;
    }
    else if (bodyName == "Earth"){
        gravitationalSelfEnergy = -4.64E-10;
    }
    else if (bodyName == "Moon"){
        gravitationalSelfEnergy = -1.88E-11;
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
        else {
            throw std::runtime_error( "Error, gravitational self energy for body "
                                      + bodyName + " can not be implemented, body radius unknown");
        }

        gravitationalSelfEnergy =
                (-3.0/5.0) *
                (gravitationalParameter*gravitationalParameter/physical_constants::GRAVITATIONAL_CONSTANT)
                /bodyRadius;
    }

    return gravitationalSelfEnergy;
}


//! Function to get the new position of the Sun due to SEP effects
//!     This function is called when creating the acceleration due to the SEP in the function below
//!     The equation is given in Genova et al. 2018, nature communication, equation (9)
//!     A simplified result has been derived and the original equations are commented out
Eigen::Vector3d getSEPCorrectedPosition(
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const NamedBodyMap& bodyMap,
        std::vector< std::string > bodyNames)
{


    // get SEP corrected position
//    Eigen::Vector3d SEPPositionCorrection;

    // Simplified solution by using taylor series expansion
    Eigen::Vector3d SEPPositionCorrectionSimplified;

    double nordtvedtParameter = ppnParameterSet->getNordtvedtParameter();

//    double nordtvedtParameter = 0.0; //for testing purposes
//    std::cout<<"eta: "<<nordtvedtParameter;

    if (nordtvedtParameter == 0.0){ // to prevent dividing by 0 in the calculations below
        SEPPositionCorrectionSimplified << 0.0, 0.0, 0.0;
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



//        double referenceCentralBodyTerm =
//                -1.0/gravitationalParameterBodyExertingAcceleration;

//        double centralBodyTerm =
//                -1.0/(
//                    gravitationalParameterBodyExertingAcceleration *
//                    (1.0 - nordtvedtParameter *
//                         gravitationalSelfEnergyBodyExertingAcceleration /
//                         ( ( gravitationalParameterBodyExertingAcceleration/physical_constants::GRAVITATIONAL_CONSTANT)
//                           * physical_constants::SPEED_OF_LIGHT*physical_constants::SPEED_OF_LIGHT
//                         )
//                    )
//                );

        // get term dependent on properties of the other bodies
//        Eigen::Vector3d summationTerm = Eigen::Vector3d::Zero();
//        Eigen::Vector3d referenceSummationTerm = Eigen::Vector3d::Zero();
        Eigen::Vector3d simplifiedSummationTerm = Eigen::Vector3d::Zero();

        const int numberOfBodies = bodyNames.size();

        for( int i = 0; i < numberOfBodies; i++ ){

           std::string currentBodyName = bodyNames.at(i);
//           Eigen::Vector3d currentBodyTerm;
//           Eigen::Vector3d currentReferenceBodyTerm;
           Eigen::Vector3d currentSimplifiedBodyTerm;

           if (currentBodyName != nameOfBodyExertingAcceleration){

               const std::shared_ptr<Body> currentBody = bodyMap.at(currentBodyName);

               Eigen::Vector3d currentBodyPosition = currentBody->getPosition();

               std::shared_ptr< gravitation::GravityFieldModel> currentBodyGravityField =
                       currentBody->getGravityFieldModel();
               double currentBodyGravitationalParameter =
                       currentBodyGravityField->getGravitationalParameter();

               double currentBodyGravitationalSelfEnergy =
                       getGravitationalSelfEnergy(currentBodyName, currentBodyGravitationalParameter);

//               currentReferenceBodyTerm =
//                       currentBodyGravitationalParameter
//                       * currentBodyPosition;

//               currentBodyTerm =
//                       (1.0 - nordtvedtParameter*
//                            currentBodyGravitationalSelfEnergy/
//                            ( ( currentBodyGravitationalParameter/physical_constants::GRAVITATIONAL_CONSTANT)
//                              * physical_constants::SPEED_OF_LIGHT*physical_constants::SPEED_OF_LIGHT
//                            )
//                       )
//                       * currentBodyGravitationalParameter
//                       * currentBodyPosition;

               currentSimplifiedBodyTerm =
                       currentBodyGravitationalParameter
                       * currentBodyPosition
                       * ( gravitationalSelfEnergyBodyExertingAcceleration
                         * 1.0/( gravitationalParameterBodyExertingAcceleration/physical_constants::GRAVITATIONAL_CONSTANT )
                         -
                         currentBodyGravitationalSelfEnergy
                         * 1.0/( currentBodyGravitationalParameter/physical_constants::GRAVITATIONAL_CONSTANT )
                       );

//               referenceSummationTerm += currentReferenceBodyTerm;
//               summationTerm += currentBodyTerm;
               simplifiedSummationTerm += currentSimplifiedBodyTerm;

           }

        }

        // multiply first and second term to get the SEP corrected position
//        SEPPositionCorrection = centralBodyTerm*summationTerm -
//                referenceCentralBodyTerm*referenceSummationTerm;

        SEPPositionCorrectionSimplified =
                (-1.0/gravitationalParameterBodyExertingAcceleration)
                * nordtvedtParameter
                * (1.0/(physical_constants::SPEED_OF_LIGHT*physical_constants::SPEED_OF_LIGHT))
                * simplifiedSummationTerm;

//        double SEPPositionCorrectionError = (SEPPositionCorrectionSimplified-SEPPositionCorrection).norm();
//        std::cout<<SEPPositionCorrectionError<< " "<<SEPPositionCorrectionError/SEPPositionCorrection.norm()<<std::endl;

    }
//    std::cout<<" dr_SEP: "<<SEPPositionCorrectionSimplified.transpose()<<std::endl;
    return SEPPositionCorrectionSimplified;
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

    Eigen::Vector3d NordtvedtPartial;
    NordtvedtPartial = Eigen::Vector3d::Zero();

    Eigen::Vector3d term1, term2, term3;
    term1 = Eigen::Vector3d::Zero();
    term2 = Eigen::Vector3d::Zero();
    term3 = Eigen::Vector3d::Zero();

    // Properties of body exerting acceleration
    Eigen::Vector3d positionBodyExertingAcceleration = bodyExertingAcceleration->getPosition( );

    std::shared_ptr< GravityFieldModel > gravityFieldCentralBody = bodyExertingAcceleration->getGravityFieldModel( );
    const double gravitationalParameterBodyExertingAcceleration =
            gravityFieldCentralBody->getGravitationalParameter();

    const double gravitationalSelfEnergyBodyExertingAcceleration =
            getGravitationalSelfEnergy(nameOfBodyExertingAcceleration,
                                       gravitationalParameterBodyExertingAcceleration);

    // Properties of body undergoing acceleration
    Eigen::Vector3d positionBodyUndergoingAcceleration = bodyUndergoingAcceleration->getPosition( );

    Eigen::Vector3d positionBodyUndergoingAccelerationWrtBodyExertingAcceleration =
            positionBodyUndergoingAcceleration - positionBodyExertingAcceleration;
    double distanceBodyUndergoingAccelerationWrtBodyExertingAcceleration =
            positionBodyUndergoingAccelerationWrtBodyExertingAcceleration.norm();

    std::shared_ptr< gravitation::GravityFieldModel> bodyUndergoingAccelerationGravityField =
            bodyUndergoingAcceleration->getGravityFieldModel();
    double gravitationalParameterBodyUndergoingAcceleration =
            bodyUndergoingAccelerationGravityField->getGravitationalParameter( );

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

        std::shared_ptr< gravitation::GravityFieldModel> currentBodyGravityField =
                currentBody->getGravityFieldModel();
        double gravitationalParameterCurrentBody =
                currentBodyGravityField->getGravitationalParameter( );

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
            term3 += gravitationalParameterCurrentBody
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

        NordtvedtPartial = term1 + term2 + term3;

//        std::cout<<NordtvedtPartial.transpose()<<" = "
//                 <<term1.transpose()<<" + "
//                 <<term2.transpose()<<" + "
//                 <<term3.transpose()<<std::endl;

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

    return NordtvedtPartial;

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


        // Set Nordtvedt parameter

//        ppnParameterSet->setNordtvedtParameter(
//                    sepViolationAccelerationSettings->nordtvedtParameter_);

        if (sepViolationAccelerationSettings->useNordtvedtConstraint_ == true){
            std::function< double( ) >  nordtvedtParameterFunction;
            nordtvedtParameterFunction =
                    std::bind( &PPNParameterSet::getNordtvedtParameterFromPpnParameters, ppnParameterSet );
        }
        else{
            std::function< double( ) >  nordtvedtParameterFunction;
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
                                sepViolationAccelerationSettings->bodyNames_);

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
                  useNordtvedtConstraintFunction
                  );


    }
    return accelerationModel;
}

}

}
