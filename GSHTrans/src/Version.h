#ifndef GSH_TRANS_VERSION_GUARD_H
#define GSH_TRANS_VERSION_GUARD_H

// The version of the library, for downstream feature tests.
//
// A dependent project may need to know which facilities are present -- the
// layered fields, the contravariant derivative and the tensor algebra all
// arrived at different points -- and a header-only library gives it no other
// way to ask. The macros here allow the question at preprocessing time:
//
//   #include <GSHTrans/src/Version.h>
//
//   #if GSHTRANS_VERSION >= GSHTRANS_VERSION_NUMBER(1, 1, 0)
//   #include <GSHTrans/Layered>
//   #endif
//
// The values are written here by hand rather than generated into the build
// tree, so that copying the GSHTrans directory into another project remains a
// complete way to vendor the library, and so that the header is the same file
// whether the library was built or merely included. CMakeLists.txt reads them
// back at configure time and refuses to configure if they disagree with the
// version declared there, so the two cannot drift apart.

#define GSHTRANS_VERSION_MAJOR 1
#define GSHTRANS_VERSION_MINOR 0
#define GSHTRANS_VERSION_PATCH 0

#define GSHTRANS_VERSION_NUMBER(major, minor, patch) \
  ((major) * 10000 + (minor) * 100 + (patch))

#define GSHTRANS_VERSION                                            \
  GSHTRANS_VERSION_NUMBER(GSHTRANS_VERSION_MAJOR,                   \
                          GSHTRANS_VERSION_MINOR,                   \
                          GSHTRANS_VERSION_PATCH)

#define GSHTRANS_VERSION_STRING "1.0.0"

namespace GSHTrans {

inline constexpr int VersionMajor = GSHTRANS_VERSION_MAJOR;
inline constexpr int VersionMinor = GSHTRANS_VERSION_MINOR;
inline constexpr int VersionPatch = GSHTRANS_VERSION_PATCH;
inline constexpr const char* VersionString = GSHTRANS_VERSION_STRING;

}  // namespace GSHTrans

#endif  // GSH_TRANS_VERSION_GUARD_H
