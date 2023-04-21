add_subdirectory(module_vdw)
add_subdirectory(module_surchem)
add_subdirectory(module_xc)
add_subdirectory(module_ewald)

list(APPEND objects
    operator.cpp
)

add_library(
    hamilt_general
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(hamilt_general)
endif()
