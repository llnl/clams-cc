cmake_minimum_required(VERSION 3.14)

# 
# CMake policy
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.30")
    cmake_policy(SET CMP0167 NEW)
endif ()

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.24")
    cmake_policy(SET CMP0135 NEW)
endif ()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------------------- #
project("Distributed-Clustering-Metrics"
	VERSION 0.1         
	DESCRIPTION "Partition comparison metrics calculated using YGM"         
	LANGUAGES CXX) 

# 
# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "CMAKE_BUILD_TYPE is set as Release")
endif()

include(FetchContent)

#
#  MPI
find_package(MPI)

# 
#  cereal 
find_package(cereal QUIET) 
if (NOT cereal_FOUND)     
	set(JUST_INSTALL_CEREAL TRUE)
	FetchContent_Declare(cereal             
		GIT_REPOSITORY https://github.com/USCiLab/cereal.git             
		GIT_TAG master             
		)     
	FetchContent_MakeAvailable(cereal) 
endif ()

# 
#  YGM
set(DESIRED_YGM_VERSION 0.9)
find_package(ygm ${DESIRED_YGM_VERSION} CONFIG)
if (NOT ygm_FOUND)
    FetchContent_Declare(
        ygm
        GIT_REPOSITORY https://github.com/LLNL/ygm
        GIT_TAG master
    )
    FetchContent_MakeAvailable(ygm)
endif ()

#
# Add executable
#
function(add_basic_executable exe_name)
	add_executable(${exe_name} ${exe_name}.cpp)
	target_include_directories(${exe_name}
                             	PRIVATE ${PROJECT_SOURCE_DIR}/include)
  	target_include_directories(${exe_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
endfunction()

#
# Link YGM
#
function(setup_ygm_target exe_name)
    target_link_libraries(${exe_name} PUBLIC MPI::MPI_CXX)
    target_link_libraries(${exe_name} PUBLIC cereal::cereal)
    target_link_libraries(${exe_name} PRIVATE ygm::ygm)
endfunction()

add_subdirectory(src)
add_subdirectory(extra)

enable_testing()
add_subdirectory(test)
