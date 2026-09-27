include_guard(GLOBAL)

function(_lmcas_architecture_array values output)
    set(_items "")
    foreach(_value IN LISTS values)
        _lmmc_quality_json_string("${_value}" _item)
        list(APPEND _items "${_item}")
    endforeach()
    list(JOIN _items "," _joined)
    set(${output} "[${_joined}]" PARENT_SCOPE)
endfunction()

function(lmcas_add_architecture_gate)
    if(NOT CMAKE_NM OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        message(FATAL_ERROR
            "Architecture gates require GCC, Clang, or AppleClang dependency files and nm")
    endif()
    if(APPLE)
        list(LENGTH CMAKE_OSX_ARCHITECTURES _lmcas_osx_architecture_count)
        if(NOT _lmcas_osx_architecture_count EQUAL 1)
            message(FATAL_ERROR
                "Darwin architecture evidence requires exactly one "
                "CMAKE_OSX_ARCHITECTURES entry: arm64 or x86_64. "
                "Disable LMCAS_ENABLE_QUALITY_GATES for universal builds.")
        endif()
        list(GET CMAKE_OSX_ARCHITECTURES 0 _lmcas_osx_architecture)
        if(NOT _lmcas_osx_architecture MATCHES "^(arm64|x86_64)$")
            message(FATAL_ERROR
                "Darwin architecture evidence supports only arm64 or x86_64; "
                "got '${_lmcas_osx_architecture}'")
        endif()
        cmake_path(GET CMAKE_NM FILENAME _lmcas_nm_name)
        if(NOT _lmcas_nm_name MATCHES "^llvm-nm(-[0-9]+(\\.[0-9]+)*)?$")
            message(FATAL_ERROR
                "Darwin architecture evidence requires llvm-nm (version-suffixed "
                "names are supported). Configure with "
                "-DCMAKE_NM=\"$(brew --prefix llvm)/bin/llvm-nm\".")
        endif()
    endif()
    set(_components "")
    set(_files "")
    foreach(_target IN LISTS LMCAS_COMPONENT_TARGETS)
        get_target_property(_rank ${_target} LMCAS_LAYER_RANK)
        get_target_property(_dependencies ${_target} LMCAS_COMPONENT_DEPENDENCIES)
        if(NOT _dependencies)
            set(_dependencies "")
        endif()
        _lmcas_architecture_array("${_dependencies}" _dependencies_json)
        list(APPEND _components
            "{\"name\":\"${_target}\",\"rank\":${_rank},\"dependencies\":${_dependencies_json}}")
        get_target_property(_sources ${_target} SOURCES)
        get_target_property(_directory ${_target} SOURCE_DIR)
        foreach(_source IN LISTS _sources)
            if(_source MATCHES "\\$<")
                message(FATAL_ERROR "Architecture ownership requires explicit files: ${_source}")
            endif()
            cmake_path(ABSOLUTE_PATH _source BASE_DIRECTORY "${_directory}" NORMALIZE)
            _lmmc_quality_json_string("${_source}" _path_json)
            list(APPEND _files "{\"path\":${_path_json},\"owner\":\"${_target}\"}")
        endforeach()
    endforeach()
    set(_tests "")
    set(_public_tests "")
    _lmmc_quality_collect_targets("${CMAKE_CURRENT_SOURCE_DIR}" _targets)
    foreach(_target IN LISTS _targets)
        get_target_property(_linkage ${_target} LMCAS_TEST_LINKAGE)
        if(NOT _linkage)
            continue()
        endif()
        get_target_property(_libraries ${_target} LINK_LIBRARIES)
        _lmcas_architecture_array("${_libraries}" _libraries_json)
        get_target_property(_test_sources ${_target} SOURCES)
        get_target_property(_test_directory ${_target} SOURCE_DIR)
        set(_absolute_sources "")
        foreach(_source IN LISTS _test_sources)
            cmake_path(ABSOLUTE_PATH _source BASE_DIRECTORY "${_test_directory}" NORMALIZE)
            list(APPEND _absolute_sources "${_source}")
        endforeach()
        _lmcas_architecture_array("${_absolute_sources}" _test_sources_json)
        if(_linkage STREQUAL "public")
            list(APPEND _public_tests "${_target}")
        endif()
        list(APPEND _tests
            "{\"name\":\"${_target}\",\"linkage\":\"${_linkage}\",\"libraries\":${_libraries_json},\"sources\":${_test_sources_json}}")
    endforeach()
    list(JOIN _components ",\n" _components_json)
    list(JOIN _files ",\n" _files_json)
    list(JOIN _tests ",\n" _tests_json)
    _lmcas_architecture_array("${CMAKE_OSX_ARCHITECTURES}" _osx_architectures_json)
    set(_properties "")
    foreach(_pair IN ITEMS
            "root|${CMAKE_CURRENT_SOURCE_DIR}"
            "build|${CMAKE_BINARY_DIR}"
            "compile_commands|${CMAKE_BINARY_DIR}/compile_commands.json"
            "generator|${CMAKE_GENERATOR}"
            "make_program|${CMAKE_MAKE_PROGRAM}"
            "nm|${CMAKE_NM}"
            "compiler|${CMAKE_CXX_COMPILER}"
            "cmake|${CMAKE_COMMAND}"
            "system_name|${CMAKE_SYSTEM_NAME}"
            "system_processor|${CMAKE_SYSTEM_PROCESSOR}")
        string(REPLACE "|" ";" _parts "${_pair}")
        list(GET _parts 0 _key)
        list(GET _parts 1 _value)
        _lmmc_quality_json_string("${_value}" _value_json)
        string(APPEND _properties "\"${_key}\":${_value_json},\n")
    endforeach()
    set(_manifest "${CMAKE_CURRENT_BINARY_DIR}/lmcas_architecture_manifest.json")
    file(WRITE "${_manifest}" "{\n${_properties}"
        "\"osx_architectures\":${_osx_architectures_json},\n"
        "\"components\":[${_components_json}],\n"
        "\"files\":[${_files_json}],\n\"tests\":[${_tests_json}]\n}\n")
    set(_command "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/check_architecture.py"
        --manifest "${_manifest}")
    add_custom_target(lmcas_architecture
        COMMAND ${_command} --report "${CMAKE_CURRENT_BINARY_DIR}/lmcas_architecture.json"
        DEPENDS ${LMCAS_COMPONENT_TARGETS} ${_public_tests}
        VERBATIM)
    if(BUILD_TESTING)
        add_test(NAME lmcas_architecture COMMAND ${_command})
        add_test(NAME lmcas_architecture_negative_fixtures
            COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/architecture_negative_fixtures.py"
                --manifest "${_manifest}")
    endif()
endfunction()
