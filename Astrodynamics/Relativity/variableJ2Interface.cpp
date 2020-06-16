/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Relativity/variableJ2Interface.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{

namespace relativity
{


// default values: correction is zero. period is solar cycle.
std::shared_ptr< VariableJ2Interface > variableJ2Interface =
        std::make_shared< VariableJ2Interface >( 0.0,
                                                 11.0*physical_constants::JULIAN_YEAR,
                                                 0.0);


}

}
