# Table of Contents
- [Table of Contents](#table-of-contents)
- [BLASW](#blasw)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

# BLASW
```BLASW``` is a header-only wrapper for ```BLAS``` and parts of the ```LAPACK``` library. It actually uses ```CBLAS``` and ```LAPACKE``` as backends. This library makes things easier by replacing the weird naming scheme of ```BLAS``` with meaningful ones and reducing the number of needed arguments for functions from 10-11 to 2-3. It also removes the need to manually find ```CBLAS``` and ```LAPACKE``` in ```CMake```.

# Installation
You can install ```BLASW``` system-wide by downloading the [latest release](https://github.com/Mathific/Blasw/releases/latest). Make sure you have ```CMake```, a ```C++``` compiler, ```CBLAS``` and ```LAPACKE``` available on your system. Then run these in repository directory:

``` shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBLASW_INSTALL=ON -DBLASW_TESTS=OFF -DCMAKE_INSTALL_PREFIX=/usr/local/ ..
sudo cmake --install .
```

The ```CMake``` options are:
+ ```BLASW_INSTALL```: Create install target for ```BLASW```.
+ ```BLASW_TESTS```: Build the tests for ```BLASW```.

Then you can use it in ```CMake```:

``` shell
find_package(BLASW REQUIRED)
add_executable(myexec main.cpp)
target_link_libraries(myexec BLASW::BLASW)
```

Or you can use ```BLASW``` as a ```CMake``` subdirectory by cloning the repository and putting it in your source directory and use it in ```CMake```:

```
add_subdirectory(Blasw/)
add_executable(myexec main.cpp)
target_link_libraries(myexec BLASW::BLASW)
```

# Usage
Documentation of ```BLASW``` API is [here](USAGE.md).

# Contributing
You can report bugs, ask questions and request features on [issues page](../../issues).

# License
This library is licensed under BSD 3-Clause permissive license. You can read it [here](LICENSE).
