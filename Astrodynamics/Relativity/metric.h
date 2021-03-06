/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_METRIC_H
#define TUDAT_METRIC_H

#include <memory>
#include <iostream>

#include <Eigen/Core>

namespace tudat
{

namespace relativity
{


//! Minkowski metric (-1,1,1,1 signature) represented as a Matrix
const static Eigen::Matrix4d minkowskiMetric = ( Eigen::Matrix4d( ) <<
                                                 -1.0, 0.0, 0.0, 0.0,
                                                 0.0, 1.0, 0.0, 0.0,
                                                 0.0, 0.0, 1.0, 0.0,
                                                 0.0, 0.0, 0.0, 1.0 ).finished( );

//! Class that stores the PPN parameters, typically used as a 'global' environment property stored in ppnParameterSet variable
class PPNParameterSet
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param parameterGamma Value of PPN parameter gamma.
     * \param parameterBeta Value of PPN parameter beta.
     */
    PPNParameterSet( const double parameterGamma,
                     const double parameterBeta,
                     const double parameterAlpha1,
                     const double parameterAlpha2,
                     const double nordtvedtParameter,
                     const double timeVaryingGravitationalParameter):
        parameterGamma_( parameterGamma ),
        parameterBeta_( parameterBeta ),
        parameterAlpha1_( parameterAlpha1 ),
        parameterAlpha2_( parameterAlpha2 ),
        nordtvedtParameter_( nordtvedtParameter ),
        timeVaryingGravitationalParameter_( timeVaryingGravitationalParameter )
    { }

    //! Destructor
    ~PPNParameterSet( ){ }

    //! Function to retrieve value of PPN parameter gamma.
    /*!
     * Function to retrieve value of PPN parameter gamma.
     * \return Value of PPN parameter gamma.
     */
    double getParameterGamma( )
    {
        return parameterGamma_;
    }

    //! Function to retrieve value of PPN parameter beta.
    /*!
     * Function to retrieve value of PPN parameter beta.
     * \return Value of PPN parameter beta.
     */
    double getParameterBeta( )
    {
        return parameterBeta_;
    }

    //! Function to retrieve value of the Nordtvedt parameter.
    double getNordtvedtParameter( )
    {
//        std::cout<<"retreived Nordtvedt parameter: "<<nordtvedtParameter_<<std::endl;
        return nordtvedtParameter_;
    }

    //! Function to retrieve value of PPN parameter alpha1.
    double getParameterAlpha1( )
    {
        return parameterAlpha1_;
    }

    //! Function to retrieve value of PPN parameter alpha1.
    double getParameterAlpha2( )
    {
        return parameterAlpha2_;
    }

    double getNordtvedtParameterFromPpnParameters( )
    {
//        std::cout<<"retreived Nordtvedt parameter from PPN parameters: "
//                 <<4.0*getParameterBeta() - getParameterGamma() - 3.0 - getParameterAlpha1() - (2.0/3.0)*getParameterAlpha2()
//                 <<std::endl;

        return 4.0*getParameterBeta() - getParameterGamma() - 3.0
                - getParameterAlpha1() - (2.0/3.0)*getParameterAlpha2();
    }

    //! Function to retrieve value of TVGP.
    double getTimeVaryingGravitationalParameter( )
    {
        return timeVaryingGravitationalParameter_;
    }

    //! Function to reset value of PPN parameter gamma.
    /*!
     * Function to reset value of PPN parameter gamma.
     * \param parameterGamma New value of PPN parameter gamma.
     */
    void setParameterGamma( const double parameterGamma )
    {
        parameterGamma_ = parameterGamma;
    }

    //! Function to reset value of PPN parameter beta.
    /*!
     * Function to reset value of PPN parameter beta.
     * \param parameterBeta New value of PPN parameter beta.
     */
    void setParameterBeta( const double parameterBeta )
    {
        parameterBeta_ = parameterBeta;
    }

    //! Function to reset value of PPN parameter alpha1.
    void setParameterAlpha1( const double parameterAlpha1 )
    {
        parameterAlpha1_ = parameterAlpha1;
    }

    //! Function to reset value of PPN parameter alpha2.
    void setParameterAlpha2( const double parameterAlpha2 )
    {
        parameterAlpha2_ = parameterAlpha2;
    }

    //! Function to reset value of Nordtvedt parameter
    void setNordtvedtParameter( const double nordtvedtParameter )
    {
        nordtvedtParameter_ = nordtvedtParameter;
    }

    //! Function to reset value of TVGP.
    void setTimeVaryingGravitationalParameter( const double timeVaryingGravitationalParameter )
    {
        timeVaryingGravitationalParameter_ = timeVaryingGravitationalParameter;
    }

protected:

    //! Value of PPN parameter gamma.
    double parameterGamma_;

    //! Value of PPN parameter beta.
    double parameterBeta_;

    //! Value of PPN parameter alpha1.
    double parameterAlpha1_;

    //! Value of PPN parameter alpha2.
    double parameterAlpha2_;

    //! Value of the Nordtvedt parameter Eta.
    double nordtvedtParameter_;

    //! Value of TimeVaryingGravitationalParameter
    double timeVaryingGravitationalParameter_;

    bool useNordtvedtConstraint_;

};

//! Global PPN parameter set, initialized upon compilation (with values equal to GR).
extern std::shared_ptr< PPNParameterSet > ppnParameterSet;

//! Global parameter denoting EP violation in proper time rate, initialized to GR value of 0 upon compilation.
extern double equivalencePrincipleLpiViolationParameter;

////! Global parameter denoting time varying gravitational parameter, initialized as zero.
//extern double timeVaryingGravitationalParameter;

}

}

#endif // TUDAT_METRIC_H
