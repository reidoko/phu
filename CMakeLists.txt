cmake_minimum_required(VERSION 3.31)

project(phu)
set(CMAKE_CXX_STANDARD 23)

set(Boost_USE_STATIC_LIBS ON)

find_package(Boost REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "Boost library not found")
endif()

add_subdirectory(src)
target_link_libraries(phu Boost::program_options)
