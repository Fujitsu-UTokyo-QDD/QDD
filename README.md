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
