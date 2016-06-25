cmake_minimum_required(VERSION 3.5)
project(ssh_craig)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp mpcomplex.cpp mpcomplex.h)
set(LIBS libgmp.so libmblas_gmp.so libmlapack_gmp.so libgmpxx.so.0 libgomp.so.1)

include_directories(/home/philipp/lib/mpack/include /usr/include)
link_directories(/home/philipp/lib/mpack/lib /usr/lib /usr/lib64)

add_executable(ssh_craig ${SOURCE_FILES})
target_link_libraries(ssh_craig ${LIBS})

find_package(Boost 1.53 COMPONENTS REQUIRED system)
find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(DOXYGEN_FOUND)