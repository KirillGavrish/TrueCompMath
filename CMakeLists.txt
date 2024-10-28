cmake_minimum_required(VERSION 3.20)

project(graphs)

if (NOT TARGET my)
    add_subdirectory(my)
endif()

file(GLOB_RECURSE GRAPH_SOURCES graphs/data_*.cpp)

foreach (GRAPH_FILE IN LISTS GRAPH_SOURCES)
    get_filename_component(GRAPH_FILE_NAME ${GRAPH_FILE} NAME_WE)
    set(EXECUTABLE_NAME make_${GRAPH_FILE_NAME})

    add_executable(${EXECUTABLE_NAME} ${GRAPH_FILE})
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE labs)
    target_include_directories(${EXECUTABLE_NAME} PUBLIC labs/src/labs)
endforeach ()
