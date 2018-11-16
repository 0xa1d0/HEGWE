cmake_minimum_required(VERSION 3.12)
project(Matcha)

find_package(CUDA QUIET)
if(CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS} .)
    SET(ALL_CUDA_LIBS ${CUDA_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_curand_LIBRARY})
    SET(LIBS ${LIBS} ${ALL_CUDA_LIBS})
    message(STATUS "CUDA_LIBRARIES: ${CUDA_INCLUDE_DIRS} ${ALL_CUDA_LIBS}")
    set(CUDA_PROPAGATE_HOST_FLAGS ON)
    set(CUDA_SEPARABLE_COMPILATION OFF)
    list( APPEND CUDA_NVCC_FLAGS -gencode=arch=compute_30,code=compute_30 )
    list( APPEND CUDA_NVCC_FLAGS -gencode=arch=compute_52,code=sm_52  )
endif()

if(CUDA_FOUND)
    #compile cuda files and add the compiled object files to your normal source files
    FILE(GLOB_RECURSE CUDA_SRC  *.cu)
    CUDA_INCLUDE_DIRECTORIES(.)
    cuda_compile(cuda_objs ${CUDA_SRC})
    SET(SRC ${SRC} ${cuda_objs})
endif()

include_directories(/home/njkolas/NTL/ntl-10.5.0/include)
include_directories(/usr/local/ssl/include)
set(LIBS ${LIBS} /home/njkolas/NTL/ntl-10.5.0/src/ntl.a)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -lgmp -m64 -march=corei7 -msse4.2 ")

set(CMAKE_CXX_STANDARD 11)

add_executable(Matcha main.cpp MtMatrix.h MtMatrix.cpp MtCiphertext.h MtCiphertext.cpp ${SRC})

target_link_libraries(Matcha ${LIBS})
