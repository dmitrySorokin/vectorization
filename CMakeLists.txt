cmake_minimum_required(VERSION 3.6.0)

project(exampleDoPri5 LANGUAGES C CXX VERSION 0.4.2)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")

set(VECCORE_INSTALL /Users/sodi/veccore_install)


find_package(VecCore REQUIRED HINTS ${VECCORE_INSTALL}/share/VecCore/cmake)
find_package(Vc REQUIRED HINTS ${VECCORE_INSTALL}/lib/cmake)

set(UMESIMD_ROOT ${VECCORE_INSTALL})
find_package(UMESIMD REQUIRED)

find_package(Geant4 REQUIRED)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cc ${sources} ${headers})

target_compile_definitions(${PROJECT_NAME} PRIVATE -DVECCORE_ENABLE_UMESIMD -DVECCORE_ENABLE_VC)
target_link_libraries(${PROJECT_NAME} PRIVATE ${Geant4_LIBRARIES} VecCore::VecCore UMESIMD::UMESIMD ${VECCORE_INSTALL}/lib/libVc.a)
