#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(RADCellSimulation)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

message("sources" ${sources} )# something changed in here for source message

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

FIND_PACKAGE( OpenMP REQUIRED) # find openMP to implement parallel computing for loops
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(radcelllib SHARED ${sources}) # creating a libarary for radcell ,changed here 
find_package(PythonLibs REQUIRED) #find python libs for calling python function in C++
include_directories(${PYTHON_INCLUDE_DIRS})

add_executable(test test.cc ${sources} ${headers})
target_link_libraries(test ${Geant4_LIBRARIES} )
target_link_libraries(test ${PYTHON_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build microdosimetry. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
file(GLOB MAC_FILES ${PROJECT_SOURCE_DIR}/*.mac 
                    ${PROJECT_SOURCE_DIR}/*.in
                    ${PROJECT_SOURCE_DIR}/*.C 
		    ${PROJECT_SOURCE_DIR}/*.py # copying the python file to the compling folder as well
                    ${PROJECT_SOURCE_DIR}/README
                    )

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
# install(TARGETS test DESTINATION bin)
install_targets(/bin test)
install_targets(/bin radcelllib )



add_subdirectory(python_wrapper) # add subdirectory for new project , changed here 
