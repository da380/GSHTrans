
add_executable(WignerExample WignerExample.cpp)
target_link_libraries(WignerExample GSHTrans)


add_executable(GridExample GridExample.cpp)
target_link_libraries(GridExample GSHTrans)

add_executable(ScalarFieldExample ScalarFieldExample.cpp)
target_link_libraries(ScalarFieldExample GSHTrans)





add_executable(Views Views.cpp)
target_link_libraries(Views GSHTrans)