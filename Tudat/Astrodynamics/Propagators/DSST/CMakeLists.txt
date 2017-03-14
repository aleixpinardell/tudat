 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #


 option(USE_NRLMSISE00 "build Tudat with NRLMSISE-00 enabled" ON)
 if(NOT USE_NRLMSISE00)
   message(STATUS "NRLMSISE-00 disabled!")
   add_definitions(-DUSE_NRLMSISE00=0)
 else()
   message(STATUS "NRLMSISE-00 enabled!")
   add_definitions(-DUSE_NRLMSISE00=1)
   # Find USE_NRLMSISE00 library on local system.
   find_package(NRLMSISE00)

   # Include NRLMSISE00 directories.
   if(NOT APPLE OR APPLE_INCLUDE_FORCE)
     include_directories(SYSTEM AFTER "${NRLMSISE00_INCLUDE_DIR}")
   else( )
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${NRLMSISE00_INCLUDE_DIR}\"")
   endif( )
 endif( )

 if(USE_CSPICE)
   list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES tudat_spice_interface cspice)
 endif()

 if(USE_NRLMSISE00)
   list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES nrlmsise00)
 endif()

 list(APPEND TUDAT_APPLICATION_PROPAGATION_LIBRARIES tudat_simulation_setup tudat_propagators
     tudat_aerodynamics tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
     tudat_electro_magnetism tudat_propulsion tudat_ephemerides tudat_numerical_integrators tudat_reference_frames
     tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_basics ${TUDAT_APPLICATION_EXTERNAL_LIBRARIES})


# Add utilities directory
# add_subdirectory( "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities" )

# Add forces directory
# add_subdirectory( "${SRCROOT}${PROPAGATORSDIR}/DSST/forces" )

# Add source files.
set(DSST_PROPAGATOR_SOURCES
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/auxiliaryElements.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/coefficientsFactories.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/osculatingToMean.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/upperBounds.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/centralBodyPertubed.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/thirdBodyPertubed.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/conservative.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/sphericHarmoicGravity.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/conservativeThirdBodyPertubed.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/nonConservative.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/zonalSphericHarmoicGravity.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/thirdBodyCentralGravity.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/atmosphericDrag.cpp"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/conservativeSolarRadiationPressure.cpp"
)

# Add header files.
set(DSST_PROPAGATOR_HEADERS
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/body.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/auxiliaryElements.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/coefficientsFactories.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/osculatingToMean.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/upperBounds.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/utilities/vectors.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/forceModel.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/centralBodyPertubed.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/thirdBodyPertubed.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/conservative.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/sphericHarmoicGravity.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/conservativeThirdBodyPertubed.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/abstract/nonConservative.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/zonalSphericHarmoicGravity.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/thirdBodyCentralGravity.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/atmosphericDrag.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/conservativeSolarRadiationPressure.h"
    "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/availableForceModels.h"
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
target_link_libraries(test_DSST ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} dsst_propagator ${Boost_LIBRARIES} )

