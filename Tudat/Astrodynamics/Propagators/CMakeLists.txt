 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Add source files.
set(PROPAGATORS_SOURCES
  "${SRCROOT}${PROPAGATORSDIR}/nBodyStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyEnckeStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/variationalEquations.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/stateTransitionMatrixInterface.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/environmentUpdateTypes.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/singleStateTypeDerivative.cpp"
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
  "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/conservativeRadiationPressure.cpp"
)

# Add header files.
set(PROPAGATORS_HEADERS
  "${SRCROOT}${PROPAGATORSDIR}/centralBodyData.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyCowellStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyEnckeStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyDSSTStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/dynamicsStateDerivativeModel.h"
  "${SRCROOT}${PROPAGATORSDIR}/singleStateTypeDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/integrateEquations.h"
  "${SRCROOT}${PROPAGATORSDIR}/bodyMassStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/variationalEquations.h"
  "${SRCROOT}${PROPAGATORSDIR}/stateTransitionMatrixInterface.h"
  "${SRCROOT}${PROPAGATORSDIR}/environmentUpdateTypes.h"
  "${SRCROOT}${PROPAGATORSDIR}/customStateDerivative.h"
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
  "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/conservativeRadiationPressure.h"
  "${SRCROOT}${PROPAGATORSDIR}/DSST/forces/availableForceModels.h"
)

# Add static libraries.
add_library(tudat_propagators STATIC ${PROPAGATORS_SOURCES} ${PROPAGATORS_HEADERS})
setup_tudat_library_target(tudat_propagators "${SRCROOT}${PROPAGATORSDIR}")

# Add unit tests.
add_executable(test_CentralBodyData "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestCentralBodyData.cpp")
setup_custom_test_program(test_CentralBodyData "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_CentralBodyData tudat_propagators tudat_basic_mathematics ${Boost_LIBRARIES})

if(USE_CSPICE)

add_executable(test_CowellStateDerivative "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestCowellStateDerivative.cpp")
setup_custom_test_program(test_CowellStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_CowellStateDerivative ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_EnckeStateDerivative "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestEnckeStateDerivative.cpp")
setup_custom_test_program(test_EnckeStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_EnckeStateDerivative ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_SequentialVariationEquationIntegration "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestSequentialVariationalEquationIntegration.cpp")
setup_custom_test_program(test_SequentialVariationEquationIntegration "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_SequentialVariationEquationIntegration ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_VariationalEquations "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestVariationalEquationPropagation.cpp")
setup_custom_test_program(test_VariationalEquations "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_VariationalEquations ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_EnvironmentModelUpdater "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestEnvironmentUpdater.cpp")
setup_custom_test_program(test_EnvironmentModelUpdater "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_EnvironmentModelUpdater ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
add_executable(test_EnvironmentModelComputations "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestEnvironmentUpdater.cpp")
setup_custom_test_program(test_EnvironmentModelComputations "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_EnvironmentModelComputations ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_BodyMassPropagation "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestBodyMassPropagation.cpp")
setup_custom_test_program(test_BodyMassPropagation "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_BodyMassPropagation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_MultiTypeStatePropagation "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestMultiTypeStatePropagation.cpp")
setup_custom_test_program(test_MultiTypeStatePropagation "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_MultiTypeStatePropagation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_DependentVariableOutput "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestDependentVariableOutput.cpp")
setup_custom_test_program(test_DependentVariableOutput "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_DependentVariableOutput ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_StoppingConditions "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestStoppingConditions.cpp")
setup_custom_test_program(test_StoppingConditions "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_StoppingConditions ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_CustomStatePropagation "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestCustomStatePropagation.cpp")
setup_custom_test_program(test_CustomStatePropagation "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_CustomStatePropagation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})


# DSST Unit Tests

add_executable(test_DSST "${SRCROOT}${PROPAGATORSDIR}/DSST/unitTests/unitTestDSST.cpp")
setup_custom_test_program(test_DSST "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_DSST ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )


endif( )
