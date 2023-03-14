pybind11_add_module(pyQDD bindings.cpp)
target_link_libraries(pyQDD PUBLIC engineSHARED)
