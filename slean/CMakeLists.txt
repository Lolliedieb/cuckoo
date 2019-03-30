cmake_minimum_required(VERSION 3.1) 
project(cuckooSlean )
set(TARGET_NAME cuckooSlean)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set (HEADERS
    sleanTrimmer.h
)

set(TRIMMER_SRC
    crypto/blake2bx.cpp	
    sleanTrimmer.cpp			
    main.cpp
)

add_executable(${TARGET_NAME} ${TRIMMER_SRC} ${HEADERS}  )

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

