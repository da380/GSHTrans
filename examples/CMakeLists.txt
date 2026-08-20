
add_executable(WignerExample WignerExample.cpp)
target_link_libraries(WignerExample GSHTrans)


add_executable(GridExample GridExample.cpp)
target_link_libraries(GridExample GSHTrans)

add_executable(FieldExample FieldExample.cpp)
target_link_libraries(FieldExample GSHTrans)

add_executable(TensorExample TensorExample.cpp)
target_link_libraries(TensorExample GSHTrans)

add_executable(wig wig.cpp)
target_link_libraries(wig GSHTrans)

add_executable(wigner3j_tests wigner3j_tests.cpp)


