#ifndef TimeVaryingGravitationalParameter_H
#define TimeVaryingGravitationalParameter_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
//#include "tudatApplications/thesis/MyApplications/TVGPInterface.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/Relativity/metric.h"


namespace tudat
{

namespace estimatable_parameters{





//! Class used to estimate PPN parameter beta
class TimeVaryingGravitationalParameter: public EstimatableParameter< double >
{

public:

    //! Constuctor
    /*!
     * Constuctor
     * \param ppnParameterSet Object used to store PPN parameters
     */
    TimeVaryingGravitationalParameter(
            const std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet
                      = relativity::ppnParameterSet  ):
        EstimatableParameter< double >( time_varying_gravitational_parameter, "global_metric" ),
      ppnParameterSet_( ppnParameterSet ){ }

    //! Destructor
    ~TimeVaryingGravitationalParameter( ) { }

    //! Function to get the current value of the PPN parameter beta.
    /*!
     * Function to get the current value of the PPN parameter beta.
     * \return Current value of the PPN parameter beta.
     */
    double getParameterValue( )
    {
        return ppnParameterSet_->getTimeVaryingGravitationalParameter();
    }

    //! Function to reset the value of the PPN parameter beta.
    /*!
     * Function to reset the value of the PPN parameter beta.
     * \param parameterValue New value of the PPN parameter beta.
     */
    void setParameterValue( double parameterValue )
    {
        ppnParameterSet_->setTimeVaryingGravitationalParameter( parameterValue );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    //! Object used to store PPN parameters
    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;

};











////! Interface class for the estimation of a gravitational parameter
//class TimeVaryingGravitationalParameter: public EstimatableParameter< double >
//{

//public:

//    //! Constructor
//    /*!
//     * Constructor
//     * \param gravityFieldModel Gravity field object containing the gravitational parameter to be estimated.
//     * \param associatedBody Name of body containing the gravityFieldModel object
//     */
//    TimeVaryingGravitationalParameter(
//            const std::shared_ptr< gravitation::TVGPModel > tvgpModel,
//            const std::string& associatedBody ):
//        EstimatableParameter< double >( time_varying_gravitational_parameter, associatedBody ),
//        tvgpModel_( tvgpModel )
//    { }

//    //! Destructor
//    ~TimeVaryingGravitationalParameter( ) { }

//    //! Function to get the current value of the gravitational parameter that is to be estimated.
//    /*!
//     * Function to get the current value of the gravitational parameter that is to be estimated.
//     * \return Current value of the gravitational parameter that is to be estimated.
//     */
//    double getParameterValue( )
//    {
//        return tvgpModel_->getTimeVaryingGravitationalParameter( );
//    }



//    //! Function to reset the value of the gravitational parameter that is to be estimated.
//    /*!
//     * Function to reset the value of the gravitational parameter that is to be estimated.
//     * \param parameterValue New value of the gravitational parameter that is to be estimated.
//     */
//    void setParameterValue( double parameterValue )
//    {
//        tvgpModel_->resetTimeVaryingGravitationalParameter( parameterValue );
//    }

//    //! Function to retrieve the size of the parameter (always 1).
//    /*!
//     *  Function to retrieve the size of the parameter (always 1).
//     *  \return Size of parameter value (always 1).
//     */
//    int getParameterSize( )
//    {
//        return 1;
//    }

//protected:

//private:

//    //! Gravity field object containing the gravitational parameter to be estimated.
//    std::shared_ptr< gravitation::TVGPModel > tvgpModel_;

//};



}

}



#endif // TimeVaryingGravitationalParameter_H
