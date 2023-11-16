# Filename Dockerfile
# This Dockerfile can be used to create an image which can be used as an environment for FoxBerry.
# The final image features:
# - OS: debian bookworm
# - libraries: gcc, cmake, ninja, boost, eigen3, OpenMPI, sqlite 3, doxygen, python, cpplint, cppcheck, Trilinos
# - Trilinos: 14.4.0 - compiled with release optimization
# - default entrypoint: /home/user
#
# Note, that this image only contains a root-user and therefore certain libraries may complain, e.g. OpenMPI. 
# To find the workarounds, search for 'execute <name of library> as root'
# To create the image, call:
#     docker build -t <name of image>:<version> .
# Once you're satisfied with the image, give it a tag, e.g.
#     docker tag <name of image> <name of repo>/<name of image>:<version>
# And finally push:
#     docker push <name of repo>/<name of image>:<version>

# A few arguments, which also can be specified by the user
ARG USER_NAME=user

# Set debian as base layer
FROM debian:12

# update system
RUN apt-get -y update
RUN apt-get -y upgrade

#install libraries
RUN apt-get install -y apt-utils
RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get install -y build-essential git cmake libboost-all-dev libopenmpi-dev libeigen3-dev libblas-dev liblapack-dev libsqlite3-dev doxygen python3 python3-pip python3-opencv cppcheck bc ninja-build rsync python-is-python3 graphviz
RUN pip install --break-system-packages cpplint virtualenv
RUN virtualenv mynotebookenv
RUN pip install --break-system-packages opencv-python jupyter jupyterlab vtk matplotlib pandas bash_kernel
RUN python -m bash_kernel.install

# create a directory for the user
WORKDIR /home/user

# clone and install google test
RUN git clone https://github.com/google/googletest.git
WORKDIR ./googletest/build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release
RUN make install -j 4

WORKDIR /home/user

# clone and install Trilinos
RUN git clone https://github.com/trilinos/Trilinos.git TrilinosGit
WORKDIR ./TrilinosGit
RUN git checkout trilinos-release-14-4-0
WORKDIR ./build
RUN cmake ..\
    -DCMAKE_CXX_COMPILER=mpic++ \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DTrilinos_USE_GNUINSTALLDIRS=TRUE \
    -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_EXTENSIONS=OFF \
    -DCMAKE_CXX_FLAGS_RELEASE_OVERRIDE="-std=c++17 -O3 -funroll-loops" \
    -DCMAKE_C_FLAGS_RELEASE_OVERRIDE="-O3 -funroll-loops" \
    -DCMAKE_Fortran_FLAGS_RELEASE_OVERRIDE="-O5 -funroll-all-loops -malign-double" \
    -DTPL_ENABLE_MPI=ON \
    -DTPL_ENABLE_gtest=OFF \
    -DTrilinos_ENABLE_OpenMP=ON \
    -DTrilinos_ENABLE_Gtest=OFF \
    -DTrilinos_ENABLE_Epetra=ON \
    -DTrilinos_ENABLE_Tpetra=ON \
    -DTrilinos_ENABLE_Xpetra=ON \
    -DTrilinos_ENABLE_AztecOO=ON \
    -DTrilinos_ENABLE_Amesos=ON \
    -DTrilinos_ENABLE_Belos=ON \
    -DTrilinos_ENABLE_Kokkos=ON \
    -DTrilinos_ENABLE_Ifpack=ON \
    -DTrilinos_ENABLE_Ifpack2=ON \
    -DTrilinos_ENABLE_Isorropia=ON \
    -DTrilinos_ENABLE_Zoltan=ON \
    -DTrilinos_ENABLE_Zoltan2=ON \
    -DTrilinos_ENABLE_Teuchos=ON \
    -DTrilinos_ENABLE_ML=ON \
    -DTrilinos_ENABLE_MueLu=ON \
    -DTpetra_ASSUME_GPU_AWARE_MPI:BOOL=0 \
    -DMueLu_ENABLE_Tutorial=OFF \
    -DTrilinos_SHOW_DEPRECATED_WARNINGS=OFF \
    -DTrilinos_HIDE_DEPRECATED_CODE=ON \
    -DCMAKE_INSTALL_PREFIX="/usr/local"

RUN make install -j 6

RUN apt-get -y update
RUN apt-get -y upgrade

# add new user
RUN useradd -m -s /bin/bash tester
USER tester
WORKDIR /home/tester
