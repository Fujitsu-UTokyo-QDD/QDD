if (BINDINGS)
  pybind11_add_module(pyQDD bindings.cpp)
  target_link_libraries(pyQDD PUBLIC engine)
endif()