add_library(AnalyzerDict SHARED)

target_include_directories(AnalyzerDict PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
)

ROOT_GENERATE_DICTIONARY(analyzer_dict Track.h SiHit.h PCHit.h CsIHit.h LINKDEF LinkDef_Analyzer.h MODULE AnalyzerDict)

target_sources(AnalyzerDict PRIVATE 
    Track.h
    Track.cpp
    SiHit.h
    PCHit.h
    CsIHit.h
)

target_link_libraries(AnalyzerDict PUBLIC ${ROOT_LIBRARIES})

set_target_properties(AnalyzerDict PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${AA_LIBRARY_DIR})

add_custom_command(TARGET AnalyzerDict POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy
	${CMAKE_CURRENT_BINARY_DIR}/libAnalyzerDict_rdict.pcm
	${AA_LIBRARY_DIR}/libAnalyzerDict_rdict.pcm
)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wvolatile HAS_VOLATILE)

if(HAS_VOLATILE)
    target_compile_options(AnalyzerDict PRIVATE -Wno-volatile)
endif()