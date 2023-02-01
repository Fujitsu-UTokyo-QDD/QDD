# QDD

## How to build & Test
### C++
```
$ cmake .-B build -DCMAKE_BUILD_TYPE=Release # Debug, RelWithDebInfo
$ cmake --build build
```

You can find executables in `build/test`.

### Python
```
$ cmake .-B build -DCMAKE_BUILD_TYPE=Release -DBINDINGS=ON
$ cmake --build build
$ cp build/qdd/py*.so . # Currently, you need to move the .so file manually.
$ poetry install
```

You can run the test as follows.
```
$ poetry run pytest
```