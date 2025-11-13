cmake_minimum_required(VERSION 3.15)

project(
	PSI
	VERSION 0.1
	DESCRIPTION "spatial partition tree"
	HOMEPAGE_URL ""
	LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # works
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --gcc-toolchain=/opt/rh/gcc-toolset-14/root/usr/")
# add_compile_options(-Wall -Wextra -Wshadow -Wconversion -Wunused -Wsign-conversion -Wfloat-equal -Wpedantic -Wuninitialized -Wformat -Wformat-security -Wmissing-braces -Wparentheses -Wsequence-point -Wswitch -Wswitch-enum -Wstrict-aliasing -Wdeprecated-declarations -Wreturn-type -Wunused-result)
# add_compile_options(-Wall -Wextra -Werror)
add_compile_options(-w)

option(DEBUG "" OFF)
option(CILKPLUS "" OFF)
option(OPENCILK "" OFF)
option(SERIAL "" OFF)
option(STDALLOC "" OFF)
option(CGAL "" OFF)
option(MEMCHECK "" OFF)
option(JEMA "" OFF)
option(CCP "" OFF)


if(CILKPLUS)
	add_compile_options(-DPARLAY_CILKPLUS -DCILK -fcilkplus)
elseif(OPENCILK)
	add_compile_options(-DPARLAY_OPENCILK -DCILK -fopencilk)
elseif(SERIAL)
	add_compile_options(-DPARLAY_SEQUENTIAL)
else()
	add_compile_options(-pthread)
endif()

if(DEBUG)
	set(CMAKE_BUILD_TYPE "Debug")
	add_compile_options(-g -march=native)
	if(NOT APPLE)
		add_compile_options(-mcx16)
	endif()
	if(CCP)
		message("CCP=ON")
		add_compile_options(-DCCP)
	endif()
else()
	set(CMAKE_BUILD_TYPE "Release")
	add_compile_options(-O3 -march=native)
	if(NOT APPLE)
		add_compile_options(-mcx16)
	endif()
endif()

if(STDALLOC)
	add_compile_options(-DPARLAY_USE_STD_ALLOC)
elseif(JEMA)
	include_directories("/usr/local/lib")
	set(CMAKE_BUILD_RPATH "/usr/local/lib")
endif()

message(STATUS "\n\n---------------------------- PARLAY -----------------------------")
add_subdirectory(include/parlaylib)
message(STATUS "---------------------------- END PARLAY -----------------------------\n\n")

message(STATUS "\n\n---------------------------- LIBMORTON -----------------------------")
set(BUILD_TESTING OFF CACHE BOOL "Disable the libmorton test")
add_subdirectory(include/libmorton)
message(STATUS "---------------------------- END LIBMORTON -----------------------------\n\n")

# NOTE: clear CMAKE_MODULE_PATH to avoid finding parlay's cmake files
list(REMOVE_ITEM CMAKE_MODULE_PATH "${PARLAY_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${PSI_SOURCE_DIR}/cmake")

# NOTE: begin of psi library
message(STATUS "\n\n---------------------------- PSI -----------------------------")

# NOTE: PSI
add_library(PSI INTERFACE)
set(PSI_INCLUDE_DIR "${PSI_SOURCE_DIR}/include/psi")
target_link_libraries(PSI INTERFACE parlay libmorton)
target_include_directories(PSI INTERFACE $<BUILD_INTERFACE:${PSI_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_subdirectory(static_analysis)
message(STATUS "---------------------------- END PSI -----------------------------\n\n")

add_executable(kd_test ${PROJECT_SOURCE_DIR}/tests/kd_test.cpp)
target_link_libraries(kd_test PUBLIC PSI)
target_include_directories(kd_test PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable(p_test ${PROJECT_SOURCE_DIR}/tests/p_test.cpp)
target_link_libraries(p_test PUBLIC PSI)
target_include_directories(p_test PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable(baselines ${PROJECT_SOURCE_DIR}/tests/baselines.cpp)
target_link_libraries(baselines PUBLIC PSI)
target_include_directories(baselines PUBLIC ${PROJECT_SOURCE_DIR}/include)

if(JEMA)
	target_link_libraries(kd_test PUBLIC libjemalloc.so)
	target_link_libraries(p_test PUBLIC libjemalloc.so)
	target_link_libraries(baselines PUBLIC libjemalloc.so)
endif()

add_executable(data_generator ${PROJECT_SOURCE_DIR}/tests/data_generate.cpp)
target_link_libraries(data_generator PUBLIC PSI)
target_include_directories(data_generator PUBLIC ${PROJECT_SOURCE_DIR}/include)


add_executable(data_washer ${PROJECT_SOURCE_DIR}/tests/data_wash.cpp)
target_link_libraries(data_washer PUBLIC PSI)
target_include_directories(data_washer PUBLIC ${PROJECT_SOURCE_DIR}/include)

find_package(Boost REQUIRED)
add_executable(boost_rtree ${PROJECT_SOURCE_DIR}/tests/boost_rtree.cpp)
target_link_libraries(boost_rtree PUBLIC PSI ${Boost_LIBRARIES})
target_include_directories(boost_rtree PUBLIC ${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS})
target_compile_options(boost_rtree PRIVATE -DBOOST_MATH_DISABLE_FLOAT128 -DBOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT)

# CGAL
if(CGAL)
	message(STATUS "\n\n---------------------------- CGAL -----------------------------")
	message(STATUS "To disable using -DCGAL=OFF")
	find_package(CGAL REQUIRED)

	# TBB
	find_package(TBB)

	# copied from CGAL default cmake
	if(TBB_FOUND AND NOT TARGET CGAL::TBB_support)

		if(NOT TARGET Threads::Threads)
			find_package(Threads REQUIRED)
		endif()

		add_library(CGAL::TBB_support INTERFACE IMPORTED)
		set_target_properties(CGAL::TBB_support PROPERTIES
			INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_TBB;NOMINMAX"
			INTERFACE_INCLUDE_DIRECTORIES "${TBB_INCLUDE_DIRS}"
			INTERFACE_LINK_LIBRARIES "TBB::tbb;TBB::tbbmalloc;Threads::Threads")
	endif()

	list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})

	# executables
	# add_executable(cgal ${PROJECT_SOURCE_DIR}/tests/cgal.cpp)
	# target_link_libraries(cgal PUBLIC PSI CGAL::CGAL CGAL::TBB_support)
	# target_include_directories(cgal PUBLIC ${PROJECT_SOURCE_DIR}/include )


	add_executable(kd_ccp ${PROJECT_SOURCE_DIR}/tests/kd_ccp.cpp)
	target_link_libraries(kd_ccp PUBLIC PSI CGAL::CGAL CGAL::TBB_support)
	target_include_directories(kd_ccp PUBLIC ${PROJECT_SOURCE_DIR}/include)

	add_executable(p_ccp ${PROJECT_SOURCE_DIR}/tests/p_ccp.cpp)
	target_link_libraries(p_ccp PUBLIC PSI CGAL::CGAL CGAL::TBB_support)
	target_include_directories(p_ccp PUBLIC ${PROJECT_SOURCE_DIR}/include)


	message(STATUS "\n\n---------------------------- END CGAL -----------------------------")
endif()

file(GLOB_RECURSE all_files
  ${PROJECT_SOURCE_DIR}/include/psi/*
)


FIND_PROGRAM(CLANG_FORMAT "clang-format")

IF(CLANG_FORMAT)
	add_custom_target(
		format
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
		COMMAND ${CLANG_FORMAT}
		-style=file
		-i
		${all_files}
	)
endif()
