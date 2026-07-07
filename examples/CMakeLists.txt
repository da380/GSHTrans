
add_executable(WignerExample WignerExample.cpp)
target_link_libraries(WignerExample GSHTrans)


add_executable(GridExample GridExample.cpp)
target_link_libraries(GridExample GSHTrans)

add_executable(FieldExample FieldExample.cpp)
target_link_libraries(FieldExample GSHTrans)

add_executable(3j 3j.cpp)
target_link_libraries(3j GSHTrans)

