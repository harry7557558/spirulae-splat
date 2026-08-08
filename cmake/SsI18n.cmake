# The translation library.
#
# A leaf: it depends on nothing and everything may depend on it. That is the
# point of giving it its own target rather than folding it into the engine --
# the backend libraries print device names and device problems, and the small
# test binaries link a backend without the engine. Without this they would
# either fail to link or have to keep their messages in English.
#
# Header-only where it can be (the catalogs are constexpr tables); this
# compiles the few things that cannot be: the current language, the locale
# chain, format() and scan().

file(GLOB SS_I18N_SOURCES CONFIGURE_DEPENDS ${SS_SRC}/i18n/*.cpp)

add_library(ss_i18n STATIC ${SS_I18N_SOURCES})
target_include_directories(ss_i18n PUBLIC ${SS_SRC})
target_compile_definitions(ss_i18n PUBLIC SS_DEFAULT_LANG=${SS_DEFAULT_LANG})
