FROM ubuntu:18.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    gnupg \
    ca-certificates \
    build-essential \
    curl \
    python-dev \
    cmake

RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/CMake-hdf5-1.10.1.tar.gz  \
    && tar -xvf CMake-hdf5-1.10.1.tar.gz  \
    && cd CMake-hdf5-1.10.1\
    && mkdir build \
    && cd build \
    && cmake ../hdf5-1.10.1 \
    &&  ../hdf5-1.10.1/configure --enable-cxx --enable-shared --prefix=/usr \
    && make -j -l6 \
    && make install

RUN curl -sL https://deb.nodesource.com/setup_10.x |  bash -
RUN apt install -y --no-install-recommends nodejs
RUN wget --no-check-certificate https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
RUN apt-get update
RUN apt-get install  -y --no-install-recommends intel-mkl-2018.2-046
RUN rm -rf /var/lib/apt/lists/*

RUN wget https://sourceforge.net/projects/arma/files/armadillo-9.200.7.tar.xz \
  && tar -xvf armadillo-9.200.7.tar.xz \
  && cd armadillo-9.200.7 \
  && cmake . \
  && ./configure \
  && make \
  && make install

WORKDIR /app
COPY . /app

RUN ["make", "build" ]
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin:/opt/intel/compilers_and_libraries/linux/lib/intel64:/usr/local/lib64:/usr/local/lib:/usr/lib:/app/LocalInstall/lib

ENTRYPOINT ["./fslam.bin"]
CMD ["100"]
