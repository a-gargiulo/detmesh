add_library(fsolve STATIC src/fsolve.c)

target_include_directories(fsolve PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# set_target_properties(fsolve PROPERTIES
    # ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_CURRENT_SOURCE_DIR}"
# )
