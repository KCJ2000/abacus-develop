add_library(
    ewald
    OBJECT
    dnrm2.cpp
    H_Ewald_pw.cpp
)

if(ENABLE_COVERAGE)
    add_coverage(ewald)
endif()

if(BUILD_TESTING)
    add_subdirectory(test)
endif()