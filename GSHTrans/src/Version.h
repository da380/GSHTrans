#ifndef GSH_TRANS_VERSION_GUARD_H
#define GSH_TRANS_VERSION_GUARD_H

/**
 * @file Version.h
 * @brief The library's version, for downstream feature tests.
 *
 * @details A dependent project may need to know which facilities are present,
 * and a header-only library gives it no other way to ask. The macros here
 * allow the question at preprocessing time:
 *
 * @code
 * #include <GSHTrans/src/Version.h>
 *
 * #if GSHTRANS_VERSION >= GSHTRANS_VERSION_NUMBER(1, 1, 0)
 * #include <GSHTrans/Layered>
 * #endif
 * @endcode
 *
 * The values are written here by hand rather than generated into the build
 * tree, so that copying the GSHTrans directory into another project remains a
 * complete way to vendor the library, and so that the header is the same file
 * whether the library was built or merely included. CMakeLists.txt reads them
 * back at configure time and refuses to configure if they disagree with the
 * version declared there, so the two cannot drift apart.
 */

/** @brief Major version. */
#define GSHTRANS_VERSION_MAJOR 1
/** @brief Minor version. */
#define GSHTRANS_VERSION_MINOR 0
/** @brief Patch version. */
#define GSHTRANS_VERSION_PATCH 0

/**
 * @brief Packs a version triple into one comparable integer.
 * @param major Major version.
 * @param minor Minor version.
 * @param patch Patch version.
 */
#define GSHTRANS_VERSION_NUMBER(major, minor, patch) \
  ((major) * 10000 + (minor) * 100 + (patch))

/** @brief This library's version, comparable against GSHTRANS_VERSION_NUMBER.
 */
#define GSHTRANS_VERSION                                                  \
  GSHTRANS_VERSION_NUMBER(GSHTRANS_VERSION_MAJOR, GSHTRANS_VERSION_MINOR, \
                          GSHTRANS_VERSION_PATCH)

/** @brief This library's version as "major.minor.patch". */
#define GSHTRANS_VERSION_STRING "1.0.0"

namespace GSHTrans {

/** @brief Major version, for use where a macro will not do. */
inline constexpr int VersionMajor = GSHTRANS_VERSION_MAJOR;
/** @brief Minor version, for use where a macro will not do. */
inline constexpr int VersionMinor = GSHTRANS_VERSION_MINOR;
/** @brief Patch version, for use where a macro will not do. */
inline constexpr int VersionPatch = GSHTRANS_VERSION_PATCH;
/** @brief Version as "major.minor.patch", for use where a macro will not do. */
inline constexpr const char* VersionString = GSHTRANS_VERSION_STRING;

}  // namespace GSHTrans

#endif  // GSH_TRANS_VERSION_GUARD_H
