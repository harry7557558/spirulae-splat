# ---------------------------------------------------------------------------
# Packaging -- the redistributable form of the app (docs/build.md#packaging).
# Never part of a normal build; ask for a target by name:
#
#     cmake --build build --target macos_app
#     cmake --build build --target macos_dmg
# ---------------------------------------------------------------------------

if(APPLE AND SS_BUILD_GUI)
    add_custom_target(macos_app
        COMMAND bash ${SS_ROOT}/tools/package_macos.sh
                --build-dir ${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${SS_ROOT}
        COMMENT "Packaging Spirula Studio.app"
        VERBATIM USES_TERMINAL)
    add_dependencies(macos_app spirula)

    add_custom_target(macos_dmg
        COMMAND bash ${SS_ROOT}/tools/package_macos.sh
                --build-dir ${CMAKE_BINARY_DIR} --dmg
        WORKING_DIRECTORY ${SS_ROOT}
        COMMENT "Packaging Spirula Studio.dmg"
        VERBATIM USES_TERMINAL)
    add_dependencies(macos_dmg spirula)
endif()
