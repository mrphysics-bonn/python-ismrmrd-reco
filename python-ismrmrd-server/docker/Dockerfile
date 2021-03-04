FROM ubuntu:20.04 AS mrd_converter

LABEL maintainer="Philipp Ehses (philipp.ehses@dzne.de)"

ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update --quiet \
    && apt-get install --no-install-recommends --no-install-suggests --yes \
       wget apt-utils apt-transport-https gnupg2 python3-dev python3-pip ca-certificates libhdf5-serial-dev h5utils cmake cmake-curses-gui libboost-all-dev libfftw3-dev git build-essential gfortran libpng-dev liblapacke-dev libpng-dev libopenblas-dev libxml2-dev libxslt1-dev 


# install mkl
COPY ./mkl.sh /tmp/
RUN /tmp/mkl.sh


RUN  mkdir -p /opt/code

RUN  cd /opt/code \
     && git clone -b bart https://github.com/pehses/python-ismrmrd-server.git


# install spiraltraj for trajectory calculation
COPY spiraltraj /opt/code/spiraltraj
RUN  cd /opt/code \
     && cd spiraltraj \
     && pip3 wheel .
     # && git clone https://github.com/mrphysics-bonn/spiral-pypulseq-example \
     # && mv spiral-pypulseq-example/spiraltraj . \


# ISMRMRD library
RUN cd /opt/code \
    && git clone https://github.com/ismrmrd/ismrmrd.git \
    && cd ismrmrd \
    && mkdir build \
    && cd build \
    && cmake ../ \
    && make -j $(nproc) \
    && make install

# siemens_to_ismrmrd converter
RUN cd /opt/code \
    && git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git \
    && cd siemens_to_ismrmrd \
    && cp /opt/code/python-ismrmrd-server/parameter_maps/* parameter_maps/ \
    && mkdir build \
    && cd build \
    && cmake ../ \
    && make -j $(nproc) \
    && make install


# compile & install the bart MRI toolbox
RUN  cd /opt/code \
     && git clone https://github.com/mrirecon/bart \
     && cd bart \
     && echo MKL=1 >> Makefile.local \
     && echo MKL_BASE=/opt/intel/mkl >> Makefile.local \
     && echo PARALLEL=1 >> Makefile.local \
     && make \
     && make install
     # && echo SLINK=1 >> Makefile.local \
     # && echo CUDA=1 >> Makefile.local \
     # && echo CUDA_BASE=/usr/local/cuda/ >> Makefile.local \
     # && echo ENABLE_MEM_CFL=1 >> Makefile.local \
     # && echo ISMRMRD=1 > Makefile.local \
     # && export ISMRMRD_HOME=/usr/local \
       

# RUN  cd /opt/code \
#      && git clone https://github.com/mrfil/PowerGrid.git \
#      && cd PowerGrid \
#      && mkdir build && cd build \
#      && apt install --no-install-recommends --no-install-suggests --yes libarmadillo-dev \
#      && cmake -DOPENACC_GPU=OFF .. \
#      && make -j $(nproc)

# ----- Start another clean build without all of the build dependencies of siemens_to_ismrmrd -----
FROM ubuntu:20.04
LABEL maintainer="Philipp Ehses (philipp.ehses@dzne.de)"

ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Copy siemens_to_ismrmrd from last stage and re-add necessary dependencies
COPY --from=mrd_converter /usr/local/lib/libismrmrd*        /usr/local/lib/
COPY --from=mrd_converter /usr/local/bin/siemens_to_ismrmrd /usr/local/bin/siemens_to_ismrmrd
COPY --from=mrd_converter /usr/local/bin/bart /usr/local/bin/bart
COPY --from=mrd_converter /opt/code/spiraltraj/SpiralTraj* /opt/code/

# Copy MKL
COPY --from=mrd_converter /opt/intel/mkl/lib /opt/intel/mkl/lib/
RUN update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so  \
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

# create pythonpath for single file packages (without requiring setup.py or __init__.py)
RUN mkdir -p /opt/code/pythonpath
ENV PYTHONPATH "${PYTHONPATH}:/opt/code/pythonpath"

# copy bart python interface to pythonpath
COPY --from=mrd_converter /opt/code/bart/python/* /opt/code/pythonpath/

# Dependencies for Python MRD server
RUN  apt-get update --quiet \
     && apt-get dist-upgrade --yes \
     && apt-get install --no-install-recommends --no-install-suggests --yes \
        ca-certificates nano git gcc python3 python3-pip h5utils libxslt1.1 libfftw3-3 liblapacke\
     && pip3 install --no-cache-dir pyxb h5py

RUN  mkdir -p /opt/code

RUN  cd /opt/code \
     && pip3 install SpiralTraj* \
     && rm SpiralTraj*

RUN  cd /opt/code \
     && git clone https://github.com/ismrmrd/ismrmrd-python.git \
     && cd /opt/code/ismrmrd-python \
     && git checkout 80fecd0 \
     && python3 setup.py install

RUN  cd /opt/code \
     && git clone https://github.com/ismrmrd/ismrmrd-python-tools.git \
     && cd /opt/code/ismrmrd-python-tools \
     && python3 setup.py install

RUN  cd /opt/code \
     && git clone -b bart https://github.com/pehses/python-ismrmrd-server.git

# Cleanup files not required after compiling
RUN  apt-get remove --yes gcc git python3-pip
RUN  apt-get autoremove --yes 
RUN  apt-get autoclean --yes
RUN  rm -r /root/.cache/pip

# add ismrmrd server start script
RUN echo "/opt/code/python-ismrmrd-server/main.py -v -H=0.0.0.0 -p=9002 -l=/tmp/share/debug/python-ismrmrd-server.log" >> /usr/bin/start_server
RUN chmod +x /usr/bin/start_server


CMD ["/usr/bin/start_server"]
