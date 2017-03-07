 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Add utilities directory
# add_subdirectory( "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities" )

# Add forces directory
# add_subdirectory( "${SRCROOT}${PROPAGATORSDIR}/DSST/forces" )

# Add source files.
set(DSST_PROPAGATOR_SOURCES
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/auxiliaryElements.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/coefficientsFactory.cpp"
#    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/jacobiPolynomials.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/osculatingToMean.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/upperBounds.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/centralBodyPertubed.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/thirdBodyPertubed.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/conservative.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/sphericHarmoicGravity.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/conservativeThirdBodyPertubed.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/zonalSphericHarmoicGravity.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/thirdBodyCentralGravity.cpp"
)

# Add header files.
set(DSST_PROPAGATOR_HEADERS
    "${SRCROOT}${PROPAGATORSDIR}/DSST/dsst.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/absoluteDate.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/auxiliaryElements.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/CjSjCoefficient.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/coefficientsFactory.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/derivativeStructure.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/DSCompiler.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/hermiteInterpolator.h"
#    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/jacobiPolynomials.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/orbit.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/osculatingToMean.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/polynomialFunction.h"
#    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/polynomialUtil.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/shortPeriodicsInterpolatedCoefficient.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/shortPeriodTerms.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/spacecraftState.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/timeSpanMap.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/upperBounds.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/forceModel.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/centralBodyPertubed.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/thirdBodyPertubed.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/conservative.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/sphericHarmoicGravity.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/conservativeThirdBodyPertubed.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/zonalSphericHarmoicGravity.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/thirdBodyCentralGravity.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/forceModelTypes.h"
)

# Add static libraries.
add_library(dsst_propagator STATIC ${DSST_PROPAGATOR_SOURCES} ${DSST_PROPAGATOR_HEADERS})
setup_tudat_library_target(dsst_propagator "${SRCROOT}${PROPAGATORSDIR}/DSST")


# Add test
macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
endmacro(setup_executable_target)

add_executable(test_DSST "${SRCROOT}${PROPAGATORSDIR}/DSST/unitTests/unitTestDSST.cpp")
setup_executable_target(test_DSST "${SRCROOT}")
target_link_libraries(test_DSST dsst_propagator tudat_propagators tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES} )

