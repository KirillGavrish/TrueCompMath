# CompMath
## There are my computational mathematics labs

### You can build this project following steps below:

```
mkdir build
cd build
cmake ..
```
### if you want to build with tests:
```
mkdir build
cd build
cmake .. -DWITH_TESTS=ON
cmake --build .
```

My code may not with old compilers, try to build this project using one of the latest ones.

For example, gcc14, by adding '-DCMAKE_CXX_COMPILER=g++-14' at building stage.

> hope it'll be useful for you
