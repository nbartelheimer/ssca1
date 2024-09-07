cmake_minimum_required(VERSION 3.5)
project(ssca1)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

option(APP_MPI "MPI")
option(APP_OSHMEM "OpenSHMEM")
option(APP_GASPI "GASPI")
option(APP_GASPI_MPI "GASPI & MPI")

include(GNUInstallDirs)

add_executable(${PROJECT_NAME} "gen_scal_data.c" "gen_sim_matrix.c" "glibc_sort.c" "main.c" "pairwise_align.c" "parameters.c" "scan_backwards.c" "util.c")

if(APP_GASPI OR APP_GASPI_MPI)
	find_package(GPI2 REQUIRED)
	if(APP_GASPI_MPI)
		find_package(MPI REQUIRED)
		set(APP_MPI "yes")
	endif()
endif()

if(APP_MPI)
	find_package(MPI REQUIRED)
endif()

if(APP_OSHMEM)
	set(CMAKE_C_COMPILER "oshcc")
	set(CMAKE_CXX_COMPILER "oshc++")
endif()

target_include_directories(${PROJECT_NAME} PRIVATE include)
target_link_libraries(${PROJECT_NAME} PUBLIC "$<$<BOOL:${APP_GASPI}>:GPI2::GPI2>" "$<$<BOOL:${APP_MPI}>:MPI::MPI_C>" "m")
target_compile_definitions(${PROJECT_NAME} PRIVATE "$<$<BOOL:${APP_GASPI}>:USE_GASPI>" "$<$<BOOL:${APP_MPI}>:USE_MPI>" "$<$<BOOL:${APP_OSHMEM}>:USE_SHMEM>")
target_compile_features(${PROJECT_NAME} PRIVATE c_std_11)
install(TARGETS ssca1 DESTINATION ${CMAKE_INSTALL_BINDIR})
