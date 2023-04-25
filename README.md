# QDD

## Build and Test
You can build the repository as follows.
```
$ cmake . -B build -DCMAKE_BUILD_TYPE=Release
$ cmake --build build -j
```
BUILD_TYPE can be either `Release`, `Debug` or `RelWithDebInfo`.

You can run the test as follows. When you make a pull request, make sure all the tests become PASSED.
```
$ ./build/test/qdd_test
```

## Python Bindings
After the above commands, try the following commands.
```
$ cp build/qdd/pyQDD.*.so .
$ poetry install
$ poetry run pytest
$ poetry build
```
You can find installable files in 'dist' directory.
