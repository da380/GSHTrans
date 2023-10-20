include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG main
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

#add_executable(TestWigner
#               TestWigner.cpp)
#target_link_libraries(TestWigner PRIVATE GSHTrans gtest_main)
include(GoogleTest)
#gtest_discover_tests(TestWigner)
