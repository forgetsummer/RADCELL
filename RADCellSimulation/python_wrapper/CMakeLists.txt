# I don't really know why these two lines are needed!
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(PYINTERFACE)

#include(${CMAKE_CURRENT_SOURCE_DIR}/FindSWIG.cmake) # this is slightly modified FindSWIG.cmake file -  had to tune it to work with OSX

FIND_PACKAGE(SWIG REQUIRED)

# INCLUDE(${SWIG_USE_FILE})
INCLUDE(UseSWIG)

FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(PythonInterp)

IF(SWIG_FOUND)
  message("FOUND SWIG")
ELSE(SWIG_FOUND)
    message("could not find SWIG in your system")
ENDIF(SWIG_FOUND)


# SET(NUMPY_INCLUDE_DIR )
# # searching for numpy include dir
# EXEC_PROGRAM(${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/core/pyinterface ARGS NumpyHeader.py
#       OUTPUT_VARIABLE NUMPY_INCLUDE_DIR)


SET(LIBS
   radcelllib   
   G4vis_management
   G4digits_hits 
   G4geometry 
   G4GMocren 
   G4materials
   G4particles 
   G4RayTracer
   G4track 
   G4visXXX
   G4error_propagation
   G4graphics_reps  
   G4modeling   
   G4persistency   
   G4readout    
   G4Tree            
   G4VRML
   G4analysis  
   G4event              
   G4gl2ps        
   G4intercoms      
   G4OpenGL     
   G4physicslists  
   G4run        
   G4visHepRep       
   G4zlib
   G4clhep     
   G4FR                 
   G4global       
   G4interfaces     
   G4parmodels  
   G4processes     
   G4tracking   
   G4vis_management
   pthread
   ${PYTHON_LIBRARIES}

)


INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include)
 
SET_SOURCE_FILES_PROPERTIES(radcell.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(radcell python radcell.i)

SWIG_LINK_LIBRARIES(radcell ${LIBS})

install(TARGETS _radcell DESTINATION bin )

set(python_files_path ${CMAKE_BINARY_DIR}/python_wrapper)

INSTALL(FILES 
    ${python_files_path}/radcell.py
	DESTINATION bin
	)



# INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
# INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
# INCLUDE_DIRECTORIES(${NUMPY_INCLUDE_DIR})
# 
# 
# SET(PYINTERFACE_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pyinterface)
# SET(CMAKE_SWIG_FLAGS "")
# 
# IF(SWIG_FOUND)
#     add_subdirectory(CC3DXML)
#     add_subdirectory(CompuCellPython)
#     
# # if(${CMAKE_SYSTEM_NAME} STREQUAL Windows)        
#     
# # else(${CMAKE_SYSTEM_NAME} STREQUAL Windows)        
#     # add_subdirectory(FiPyInterface)
# # endif(${CMAKE_SYSTEM_NAME} STREQUAL Windows)            
# 
# #    add_subdirectory(PlayerPythonOld)
# 
#     FIND_PACKAGE(VTK REQUIRED)
#     message("VTK_MAJOR_VERSION=" ${VTK_MAJOR_VERSION})
#     if (${VTK_MAJOR_VERSION} EQUAL "6")    
#         add_definitions(-DVTK6)
#     endif()    
#     
# 
#     add_subdirectory(PlayerPythonNew)
#     add_subdirectory(SerializerDE)
#     add_subdirectory(PyPlugin)
#     add_subdirectory(Example)
#     add_subdirectory(Fields)
# if(USE_DOLFIN)    
#     add_subdirectory(dolfinCC3D)
# endif(USE_DOLFIN)        
    

#CMAKE_MINIMUM_REQUIRED(VERSION 2.2)
#PROJECT(PYINTERFACE)
