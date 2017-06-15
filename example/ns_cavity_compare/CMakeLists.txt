cmake_minimum_required(VERSION 2.8)
project(carpio_1d_poisson)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/)

set(SOURCE_FILES main.cpp)

set(PROJECT_SOURCE_DIR ./)

set(DIR_CARPIO ../..)


# location of header files
include_directories(
    ${DIR_CARPIO}/lib
)
message(STATUS "INCLUDE_DIR=${DIR_CARPIO}")

set(Python_ADDITIONAL_VERSIONS 3.4)
#add_subdirectory(${PROJECT_SOURCE_DIR}/lib/ts)
#add_subdirectory(${DIR_CARPIO}/lib/utility)
#add_subdirectory(${DIR_CARPIO}/lib/io)

# python =======
set(REQUIRED_PYTHON_VERSION 3)
message(STATUS "======== LOOKING FOR ${REQUIRED_PYTHON_VERSION} ========================")
find_package ( PythonInterp ${REQUIRED_PYTHON_VERSION} REQUIRED)
find_package ( PythonLibs ${REQUIRED_PYTHON_VERSION} REQUIRED)
#find_package ( PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED )


message(STATUS "PYTHON_EXECUTABLE:FILEPATH=${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_LIBRARY:FILEPATH=${PYTHON_LIBRARY}")
message(STATUS "PYTHON_LIBRARIES:=${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_INCLUDE_DIR:FILEPATH=${PYTHON_INCLUDE_DIR}")
#message(STATUS "PYTHON_FRAMEWORK_INCLUDES=${PYTHON_FRAMEWORK_INCLUDES}")
message(STATUS "PYTHONLIBS_VERSION_STRING=${PYTHONLIBS_VERSION_STRING}")
#message(STATUS "PYTHON_FRAMEWORKS=${Python_FRAMEWORKS}")

IF(PYTHONLIBS_FOUND)
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_DIRS}")
ELSE()
  MESSAGE(FATAL_ERROR "Unable to find PythonLibs.")
ENDIF()


#unset(PYTHON_EXECUTABLE CACHE)
#unset(PYTHON_LIBRARY CACHE)
#unset(PYTHON_INCLUDE_DIR CACHE)
# end python ===



find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(main ${SOURCE_FILES})


add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/utility ${CMAKE_CURRENT_BINARY_DIR}/lib/utility)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/io ${CMAKE_CURRENT_BINARY_DIR}/lib/io)


# link the libraries to the executable
target_link_libraries (main io)
target_link_libraries (main utility)
target_link_libraries (main gtest)
target_link_libraries (main glog)
target_link_libraries (main pthread)

if(PYTHON_LIBRARIES)
  TARGET_LINK_LIBRARIES(main ${PYTHON_LIBRARIES})
ELSE()
  MESSAGE(FATAL_ERROR "Unable to link python libraries")
ENDIF()

if(VTK_LIBRARIES)
    target_link_libraries(main ${VTK_LIBRARIES})
else()
    target_link_libraries(main vtkHybrid vtkWidgets)
endif()
