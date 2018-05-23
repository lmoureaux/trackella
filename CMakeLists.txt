cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# Use our own FindROOT, not the system one
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Require C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Epiphany
set(ESDK $ENV{EPIPHANY_HOME})
set(ELIBS ${ESDK}/tools/host/lib)
set(EINCS ${ESDK}/tools/host/include)
link_directories(find_doublets PRIVATE ${ELIBS})

# Dependencies
find_package(ROOT 6 REQUIRED COMPONENTS Table TreePlayer)

add_executable(print_event_stats
    src/print_event_stats.cpp
    src/doublet_finder.cpp
    src/eventreader.cpp
)
target_include_directories(print_event_stats SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(print_event_stats PUBLIC ${ROOT_LIBRARIES})

add_executable(find_doublets
    src/epiphany_doublet_finder.cpp
    src/find_doublets.cpp
    src/doublet_finder.cpp
    src/eventreader.cpp
)
target_include_directories(find_doublets SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(find_doublets PUBLIC ${ROOT_LIBRARIES})

target_include_directories(find_doublets SYSTEM PRIVATE ${EINCS})
target_link_libraries(find_doublets PRIVATE e-hal e-loader)

add_executable(asm-test src/asm-test.cpp)
