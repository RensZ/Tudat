
#ifndef TUDAT_STRONGEQUIVALENCEPRINCIPLE_H
#define TUDAT_STRONGEQUIVALENCEPRINCIPLE_H

#include <Tudat/SimulationSetup/tudatEstimationHeader.h>

namespace tudat
{

namespace relativity
{



//! get the gravitational self energy of a body, which is needed in the function below
//! values and equation come from Genova et al. 2018, nature communications
//! radii from //https://nssdc.gsfc.nasa.gov/planetary/factsheet/<INSERTPLANET>fact.html
double getGravitationalSelfEnergy(
        const std::string bodyName,
        const double gravitationalParameter);


//! Function to get the new position of the Sun due to SEP effects
//!     This function is called when creating the acceleration due to the SEP in the function below
//!     The equation is given in Genova et al. 2018, nature communication, equation (9)
//!     A simplified result has been derived and the original equations are commented out
Eigen::Vector3d getSEPCorrectedPosition(
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const NamedBodyMap& bodyMap,
        std::vector< std::string > bodyNames);



//! the partial equation of the SEP violation acceleration w.r.t. the Nordtvedt parameter is calculated
//!     This is passed on via sepViolationAcceleration to sepViolationAccelerationPartial
//!     And is calculated using Genova et al. 2018, nature communication, equation (10)
Eigen::Vector3d getNordtvedtPartial(
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const NamedBodyMap& bodyMap,
        std::vector< std::string > bodyNames);

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEACCELERATIONMODELS_H

