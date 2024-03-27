rm -f CMakeCache.txt cmake_install.cmake Makefile CTestTestfile.cmake
rm -rf CMakeFiles
rm -rf _deps
find src qdd test -type f \( -name 'CMakeCache.txt' -o -name 'cmake_install.cmake' -o -name 'Makefile' \) -exec rm {} +
find src qdd test -type d -name 'CMakeFiles' -exec rm -r {} +
git clean -ffdx . -e dist -e qdd.egg-info
