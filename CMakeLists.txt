cmake_minimum_required(VERSION 3.15)

project(
	PSPT
	VERSION 0.1
	DESCRIPTION "spatial partition tree"
	HOMEPAGE_URL ""
	LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # works
# add_compile_options(-Wall -Wextra -Wshadow -Wconversion -Wunused -Wsign-conversion -Wfloat-equal -Wpedantic -Wuninitialized -Wformat -Wformat-security -Wmissing-braces -Wparentheses -Wsequence-point -Wswitch -Wswitch-enum -Wstrict-aliasing -Wdeprecated-declarations -Wreturn-type -Wunused-result)
add_compile_options(-Wall -Wextra -Werror)

option(DEBUG "" OFF)
option(CILKPLUS "" OFF)
option(OPENCILK "" OFF)
option(SERIAL "" OFF)
option(STDALLOC "" OFF)
option(CGAL "" ON)
option(MEMCHECK "" OFF)
option(JEMA "" ON)
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
	add_compile_options(-g -mcx16 -march=native)
	if(CCP)
		message("CCP=ON")
		add_compile_options(-DCCP)
	endif()
else()
	set(CMAKE_BUILD_TYPE "Release")
	add_compile_options(-O3 -mcx16 -march=native)
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

# NOTE: clear CMAKE_MODULE_PATH to avoid finding parlay's cmake files
list(REMOVE_ITEM CMAKE_MODULE_PATH "${PARLAY_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${PSPT_SOURCE_DIR}/cmake")

add_library(PSPT INTERFACE)
set(PSPT_INCLUDE_DIR "${PSPT_SOURCE_DIR}/include/pspt")
target_link_libraries(PSPT INTERFACE parlay)
target_include_directories(PSPT INTERFACE $<BUILD_INTERFACE:${PSPT_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_subdirectory(static_analysis)

include(GNUInstallDirs)
install(TARGETS PSPT
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_executable(test ${PROJECT_SOURCE_DIR}/tests/test.cpp)
target_link_libraries(test PUBLIC PSPT)
target_include_directories(test PUBLIC ${PROJECT_SOURCE_DIR}/include)

if(JEMA)
	target_link_libraries(test PUBLIC libjemalloc.so)
endif()

add_executable(data_generator ${PROJECT_SOURCE_DIR}/tests/data_generate.cpp)
target_link_libraries(data_generator PUBLIC PSPT)
target_include_directories(data_generator PUBLIC ${PROJECT_SOURCE_DIR}/include)

find_package(Boost REQUIRED)
add_executable(rtree ${PROJECT_SOURCE_DIR}/tests/rtree.cpp)
target_link_libraries(rtree PUBLIC PSPT ${Boost_LIBRARIES})
target_include_directories(rtree PUBLIC ${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS})

# CGAL
if(CGAL)
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
	# target_link_libraries(cgal PUBLIC PSPT CGAL::CGAL CGAL::TBB_support)
	# target_include_directories(cgal PUBLIC ${PROJECT_SOURCE_DIR}/include )

	add_executable(ccp ${PROJECT_SOURCE_DIR}/tests/ccp.cpp)
	target_link_libraries(ccp PUBLIC PSPT CGAL::CGAL CGAL::TBB_support)
	target_include_directories(ccp PUBLIC ${PROJECT_SOURCE_DIR}/include)
endif()

set(all_files
  ${PROJECT_SOURCE_DIR}/include/pspt/*
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
