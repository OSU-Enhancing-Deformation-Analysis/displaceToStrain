# CMakeList.txt : CMake project for strain_calc, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("strain_calc")

# set opencv directory path
# set(OpenCV_DIR "C:/Users/Alex/Downloads/opencv/build")

#find opencv package
find_package( OpenCV REQUIRED)

# Add source to this project's executable.
add_executable (strain_calc "strain_calc.cpp" "strain_calc.h" "cnpy.h" "cnpy.cpp")

target_include_directories(strain_calc 
    PRIVATE 
    ${CMAKE_SOURCE_DIR}/dependencies/ 
    ${OpenCV_INCLUDE_DIRS}
)

target_link_libraries( strain_calc PRIVATE ${OpenCV_LIBS} )

if (CMAKE_VERSION VERSION_GREATER 3.12)
set_property(TARGET strain_calc PROPERTY CXX_STANDARD 20)
endif()