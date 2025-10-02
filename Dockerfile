# Filename Dockerfile
# This Dockerfile can be used to create an image which can be used as an environment for FoxBerry.
# The final image features:
# - OS: debian bookworm
# - libraries: gcc, cmake, clang, ninja, boost, eigen3, OpenMPI, sqlite 3, doxygen, python, cpplint, cppcheck, Trilinos
# - Trilinos: 16.1.0 - compiled with release optimization
# - default entrypoint: /home/tester
#
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

# Add non-free reposity source for CUDA
RUN echo 'deb http://deb.debian.org/debian bookworm main contrib non-free non-free-firmware' >> /etc/apt/sources.list

# update system
RUN apt-get -y update
RUN apt-get -y upgrade

#install libraries
RUN apt-get install -y apt-utils
RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get install -y build-essential git cmake wget software-properties-common libboost-all-dev libopenmpi-dev libeigen3-dev libblas-dev liblapack-dev libsqlite3-dev libgtest-dev bc ninja-build rsync graphviz mesa-common-dev mesa-utils freeglut3-dev
RUN apt-get install -y clang clang-format clangd clang-tidy nvidia-cuda-dev nvidia-cuda-toolkit libomp-dev
RUN apt-get install -y doxygen cpplint python3 python3-pip python3-opencv cppcheck python-is-python3
RUN pip install --break-system-packages cppcheck-junit cpplint-junit doxygen-junit numpy pandas matplotlib vtk compdb clang-tidy

# create a directory for the user
WORKDIR /home/user

# Silence git safe.directory warnings
RUN git config --add --system safe.directory '*'

WORKDIR /home/user

# clone and install Trilinos with gcc
RUN git clone https://github.com/trilinos/Trilinos.git TrilinosGit
WORKDIR ./TrilinosGit
RUN git pull
RUN git checkout trilinos-release-16-1-0
WORKDIR ./build_gcc
RUN cmake ..\
    -GNinja \
    -DCMAKE_CXX_COMPILER=mpic++ \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_CXX_STANDARD=20 \
    -DTrilinos_USE_GNUINSTALLDIRS=TRUE \
    -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_EXTENSIONS=OFF \
    -DTPL_ENABLE_MPI=ON \
    -DTPL_ENABLE_gtest=OFF \
    -DTrilinos_ENABLE_OpenMP=ON \
    -DTrilinos_ENABLE_Gtest=OFF \
    -DTrilinos_ENABLE_Tpetra=ON \
    -DTrilinos_ENABLE_Xpetra=ON \
    -DTrilinos_ENABLE_Amesos2=ON \
    -DTrilinos_ENABLE_Belos=ON \
    -DTrilinos_ENABLE_Kokkos=ON \
    -DTrilinos_ENABLE_Ifpack2=ON \
    -DTrilinos_ENABLE_Zoltan=ON \
    -DTrilinos_ENABLE_Zoltan2=ON \
    -DTrilinos_ENABLE_Teuchos=ON \
    -DTrilinos_ENABLE_MueLu=ON \
    -DTpetra_ASSUME_GPU_AWARE_MPI:BOOL=0 \
    -DMueLu_ENABLE_Tutorial=OFF \
    -DTrilinos_SHOW_DEPRECATED_WARNINGS=OFF \
    -DTrilinos_HIDE_DEPRECATED_CODE=ON \
    -DCMAKE_INSTALL_PREFIX="/usr/local/trilinos_gcc"

RUN ninja install -j 6

# building trilinos with clang
WORKDIR /home/user/TrilinosGit/build_clang
ENV CXX=mpic++
ENV OMPI_CC=clang
ENV OMPI_CXX=clang++
RUN cmake ..\
    -GNinja \
    -DCMAKE_CXX_COMPILER=mpic++ \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_CXX_STANDARD=20 \
    -DTrilinos_USE_GNUINSTALLDIRS=TRUE \
    -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_EXTENSIONS=OFF \
    -DTPL_ENABLE_MPI=ON \
    -DTPL_ENABLE_gtest=OFF \
    -DTrilinos_ENABLE_OpenMP=ON \
    -DTrilinos_ENABLE_Gtest=OFF \
    -DTrilinos_ENABLE_Tpetra=ON \
    -DTrilinos_ENABLE_Xpetra=ON \
    -DTrilinos_ENABLE_Amesos2=ON \
    -DTrilinos_ENABLE_Belos=ON \
    -DTrilinos_ENABLE_Kokkos=ON \
    -DTrilinos_ENABLE_Ifpack2=ON \
    -DTrilinos_ENABLE_Zoltan=ON \
    -DTrilinos_ENABLE_Zoltan2=ON \
    -DTrilinos_ENABLE_Teuchos=ON \
    -DTrilinos_ENABLE_MueLu=ON \
    -DTpetra_ASSUME_GPU_AWARE_MPI:BOOL=0 \
    -DMueLu_ENABLE_Tutorial=OFF \
    -DTrilinos_SHOW_DEPRECATED_WARNINGS=OFF \
    -DTrilinos_HIDE_DEPRECATED_CODE=ON \
    -DCMAKE_INSTALL_PREFIX="/usr/local/trilinos_clang"

RUN ninja install -j 6

#resetting environment variables
ENV OMPI_CC=gcc
ENV OMPI_CXX=g++

# install vtk
WORKDIR /home/user
RUN git clone https://gitlab.kitware.com/vtk/vtk.git
WORKDIR ./vtk
RUN git checkout v9.3.0
WORKDIR ./build
RUN cmake ..\
     -GNinja \
     -DCMAKE_BUILD_TYPE=Release\
     -DVTK_USE_MPI=ON\
     -DVTK_USE_CUDA=OFF\
     -DVTK_SMP_IMPLEMENTATION_TYPE=OpenMP\
     -DVTK_LEGACY_REMOVE=ON\
     -DVTK_USE_FUTURE_CONST=ON\
     -DVTK_USE_FUTURE_BOOL=ON

RUN ninja install -j 6

RUN apt-get -y update
RUN apt-get -y upgrade

# add new user
RUN useradd -m -s /bin/bash tester
USER tester
WORKDIR /home/tester
