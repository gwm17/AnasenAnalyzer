add_subdirectory(Dict)

add_executable(AnasenAnalyzer)

target_include_directories(AnasenAnalyzer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../vendor/catima/
    ${CMAKE_CURRENT_SOURCE_DIR}/../vendor/fmt/include/
    ${ROOT_INCLUDE_DIRS}
)

target_sources(AnasenAnalyzer PRIVATE
    Analyzer/MassLookup.h
    Analyzer/MassLookup.cpp
    Analyzer/Target.h
    Analyzer/Target.cpp
    Analyzer/CutHandler.h
    Analyzer/CutHandler.cpp
    Analyzer/GateHandler.h
    Analyzer/Reconstructor.h
    Analyzer/Reconstructor.cpp
    Analyzer/Analyzer.h
    Analyzer/Analyzer.cpp
    Analyzer/Application.h
    Analyzer/Application.cpp
    main.cpp
)

target_link_libraries(AnasenAnalyzer PRIVATE AnalyzerDict ${ROOT_LIBRARIES} fmt::fmt catima)

set_target_properties(AnasenAnalyzer PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${AA_BINARY_DIR})

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wvolatile HAS_VOLATILE)

if(HAS_VOLATILE)
    target_compile_options(AnasenAnalyzer PRIVATE -Wno-volatile)
endif()