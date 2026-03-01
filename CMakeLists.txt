cmake_minimum_required(VERSION 3.16)
project(hex-planets C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(SOKOL_DIR ${CMAKE_SOURCE_DIR}/third_party/sokol)

# Platform detection for Sokol backend
if(EMSCRIPTEN)
    set(SOKOL_BACKEND "SOKOL_GLES3")
elseif(WIN32)
    set(SOKOL_BACKEND "SOKOL_D3D11")
elseif(APPLE)
    set(SOKOL_BACKEND "SOKOL_METAL")
else()
    set(SOKOL_BACKEND "SOKOL_GLCORE")
endif()

# Sokol implementation library (single compilation unit)
add_library(sokol STATIC src/sokol_impl.c)
target_include_directories(sokol PUBLIC ${SOKOL_DIR})
target_compile_definitions(sokol PUBLIC ${SOKOL_BACKEND})

if(WIN32)
    target_link_libraries(sokol PUBLIC kernel32 user32 gdi32 d3d11 dxgi dxguid ole32 dbghelp)
elseif(APPLE)
    target_link_libraries(sokol PUBLIC
        "-framework Metal" "-framework Cocoa" "-framework QuartzCore"
        "-framework MetalKit")
elseif(UNIX AND NOT EMSCRIPTEN)
    target_link_libraries(sokol PUBLIC X11 Xi Xcursor GL dl pthread m)
endif()

# Main executable
add_executable(hex-planets
    src/main.c
    src/planet.c
    src/mesh.c
    src/camera.c
    src/render.c
    src/atmosphere.c
    src/lod.c
    src/hex_terrain.c
    src/job_system.c
    src/math_utils.c
    src/screenshot.c
    src/crash_handler.c
    src/config.c
)

target_include_directories(hex-planets PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/third_party
    ${CMAKE_SOURCE_DIR}/shaders
    ${SOKOL_DIR}
)

target_link_libraries(hex-planets PRIVATE sokol)

if(WIN32)
    set_target_properties(hex-planets PROPERTIES WIN32_EXECUTABLE TRUE)
endif()

# ---- Shader compilation with sokol-shdc ----
set(SOKOL_SHDC "${CMAKE_SOURCE_DIR}/tools/sokol-shdc${CMAKE_EXECUTABLE_SUFFIX}")
set(SHADER_SLANG "hlsl5:glsl430:glsl300es:wgsl")
set(SHADER_DIR "${CMAKE_SOURCE_DIR}/shaders")

file(GLOB SHADER_SOURCES "${SHADER_DIR}/*.glsl")
set(SHADER_HEADERS)

foreach(SHADER_SRC ${SHADER_SOURCES})
    get_filename_component(SHADER_NAME ${SHADER_SRC} NAME)
    set(SHADER_OUT "${SHADER_DIR}/${SHADER_NAME}.h")
    list(APPEND SHADER_HEADERS ${SHADER_OUT})

    add_custom_command(
        OUTPUT ${SHADER_OUT}
        COMMAND ${SOKOL_SHDC}
            --input ${SHADER_SRC}
            --output ${SHADER_OUT}
            --slang ${SHADER_SLANG}
        DEPENDS ${SHADER_SRC}
        COMMENT "Compiling shader: ${SHADER_NAME}"
        VERBATIM
    )
endforeach()

add_custom_target(shaders DEPENDS ${SHADER_HEADERS})
add_dependencies(hex-planets shaders)

# Emscripten web build settings
if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
    target_link_options(hex-planets PRIVATE
        -sUSE_WEBGL2=1
        -sALLOW_MEMORY_GROWTH=1
        -sINITIAL_MEMORY=67108864
        -sTOTAL_STACK=8388608
        --shell-file=${CMAKE_SOURCE_DIR}/web/shell.html
    )
    # Copy output to dist/ after build
    add_custom_command(TARGET hex-planets POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/dist
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:hex-planets> ${CMAKE_SOURCE_DIR}/dist/index.html
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/hex-planets.js ${CMAKE_SOURCE_DIR}/dist/hex-planets.js
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/hex-planets.wasm ${CMAKE_SOURCE_DIR}/dist/hex-planets.wasm
        COMMENT "Copying web build to dist/"
    )
endif()

if(MSVC)
    target_compile_options(hex-planets PRIVATE /W3)
    target_compile_options(sokol PRIVATE /W3)
else()
    target_compile_options(hex-planets PRIVATE -Wall -Wextra)
    target_compile_options(sokol PRIVATE -Wall -Wextra)
endif()

# ---- Test executable (no GPU required) ----
if(NOT EMSCRIPTEN)
    add_executable(hex-tests
        tests/test_main.c
        tests/test_math.c
        tests/test_planet.c
        tests/test_mesh.c
        src/planet.c
        src/mesh.c
        src/math_utils.c
    )
    target_include_directories(hex-tests PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/third_party
        ${CMAKE_SOURCE_DIR}/tests
    )
    if(MSVC)
        target_compile_options(hex-tests PRIVATE /W3)
    else()
        target_compile_options(hex-tests PRIVATE -Wall -Wextra)
        target_link_libraries(hex-tests PRIVATE m)
    endif()
endif()
