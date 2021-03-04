#!/bin/bash

# install mkl
apt update && apt install -y --force-yes apt-transport-https && \
  wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB && \
  apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB && \
  sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list' && \
  apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install cpio intel-mkl-64bit-2018.3-051 && \
  (find /opt/intel -name "ia32*" -exec rm -rf {} \; || echo "removing ia32 binaries") ; \
  (find /opt/intel -name "examples" -type d -exec rm -rf {} \; || echo "removing examples") ; \
  (find /opt/intel -name "benchmarks" -exec rm -rf {} \; || echo "removing benchmarks") ; \
  (find /opt/intel -name "documentation*" -exec rm -rf {} \; || echo "removing documentation") ; \
  (rm -rf /opt/intel/mkl/interfaces ) ; \
  (rm -rf /opt/intel/mkl/lib/intel64/*.a ) ; \
  (rm -rf /opt/intel/mkl/lib/intel64/*mpi*.so ) ; \
  (rm -rf /opt/intel/mkl/lib/intel64/*tbb*.so ) ; \
  (rm -rf /opt/intel/mkl/lib/intel64/*pgi*.so ) ; \
  (rm -rf /opt/intel/mkl/lib/intel64/*mc*.so ) ; \
  (rm -rf /opt/intel/mkl/lib/intel64/*blacs*.so ) ; \
  (rm -rf /opt/intel/mkl/lib/intel64/*scalapack*.so ) ; \
  (rm -rf /opt/intel/mkl/lib/intel64/*gf*.so ) ; \
  (rm -rf /opt/intel/mkl/lib/intel64/*mic*.so ) ; \
  # apt purge intel-tbb* intel-psxe* && \
  apt-get clean autoclean && \
  apt-get autoremove -y && \
  ln -s -f bash /bin/sh && \
  rm -rf /usr/share/doc && \
  echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/intel.conf && \
  ldconfig && \
  echo "source /opt/intel/mkl/bin/mklvars.sh intel64" >> /etc/bash.bashrc

update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so  \
    libblas.so-x86_64-linux-gnu      /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 && \
  update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so.3  \
    libblas.so.3-x86_64-linux-gnu    /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 && \
  update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so   \
    liblapack.so-x86_64-linux-gnu    /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 && \
  update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so.3 \
    liblapack.so.3-x86_64-linux-gnu  /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 && \
  echo "/opt/intel/lib/intel64"     >  /etc/ld.so.conf.d/mkl.conf && \
  echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/mkl.conf && \
  ldconfig && \
  echo "MKL_THREADING_LAYER=GNU" >> /etc/environment
