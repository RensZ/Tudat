/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_VARIABLEJ2INTERFACE_H
#define TUDAT_VARIABLEJ2INTERFACE_H

#include <memory>
#include <iostream>
#include <iomanip>

#include <Eigen/Core>

namespace tudat
{

namespace relativity
{


//! Class that stores the PPN parameters, typically used as a 'global' environment property stored in ppnParameterSet variable
class VariableJ2Interface
{
public:

    VariableJ2Interface( const double amplitude,
                         const double period,
                         const double phase,
                         const double meanJ2):
        amplitude_( amplitude ),
        period_( period ),
        phase_( phase ),
        meanJ2_( meanJ2 )
    { }

    //! Destructor
    ~VariableJ2Interface( ){ }

    double getAmplitude( ){
        std::cout<<std::setprecision(15)<<"J2 amplitude is retreived: "<<amplitude_<<std::setprecision(5)<<std::endl;
        return amplitude_; }
    double getPeriod( ){return period_; }
    double getPhase( ){return phase_; }
    double getMeanJ2( ){
        //std::cout<<std::setprecision(15)<<"J2 mean is retreived: "<<meanJ2_<<std::setprecision(5)<<std::endl;
        return meanJ2_; }

    void setAmplitude( const double amplitude ){
        std::cout<<std::setprecision(15)<<"J2 amplitude gets an update. old value: "<<amplitude_<<" new value: "<<amplitude<<std::setprecision(5)<<std::endl;
        amplitude_ = amplitude; }
    void setPeriod( const double period ){ period_ = period; }
    void setPhase( const double phase ){ phase_ = phase; }
    void setMeanJ2( const double meanJ2 ){
        std::cout<<std::setprecision(15)<<"J2 mean gets an update. old value: "<<meanJ2_<<" new value: "<<meanJ2<<std::setprecision(5)<<std::endl;
        meanJ2_ = meanJ2; }

protected:

    double amplitude_;
    double period_;
    double phase_;
    double meanJ2_;

};

//! Global PPN parameter set, initialized upon compilation (with values equal to GR).
extern std::shared_ptr< VariableJ2Interface > variableJ2Interface;

}

}

#endif // TUDAT_VARIABLEJ2INTERFACE_H
