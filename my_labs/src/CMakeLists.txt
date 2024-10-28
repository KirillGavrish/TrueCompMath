file(GLOB_RECURSE SOURCES labs/*.h labs/*.c labs/*.hpp labs/*.cpp)
add_library(labs INTERFACE ${SOURCES})

target_link_libraries(labs INTERFACE eigen)

target_include_directories(labs INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(labs INTERFACE cxx_std_20)
