# STITCH
Linking regulatory regions to genes

## Building the tool
To build the **STICH** *cmake*, a C++11 compiler, and the boost library must be available.

Make sure you are in the projects main directory. Then, generate a build folder by typing:
`mkdir build`

Change into that folder:
`cd build`

To run *cmake* type:
`cmake ..`

To finally build the project type:
`make`

You can speed up the build process by using the *-j* option, specifying how many cores should be used for building

To execute all tests type:
`make test`

Tests can also be executed individually. They are located in the folder
`build/test`
