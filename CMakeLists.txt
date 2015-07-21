project(eLISA_DPC)
cmake_minimum_required(VERSION 2.6)
#SET (CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS_INIT} $ENV{LDFLAGS}
     CACHE STRING "Flags used by the linker during the creation of dll's.")
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules)
MESSAGE(STATUS “${CMAKE_MODULE_PATH} READING CONFIGURATION”)
add_subdirectory(SIM)