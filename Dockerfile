FROM ubuntu:20.04

# this is for timezone config
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Amsterdam
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install packages
RUN apt-get update
RUN apt-get install -y apt-utils
RUN apt-get install -y build-essential
RUN apt-get install -y libboost-all-dev
RUN apt-get install -y liblapack-dev libblas-dev
RUN apt-get install -y automake
RUN apt-get install -y cmake
RUN apt-get install -y make 
RUN apt-get install -y g++-10 
RUN apt-get install -y libgmp-dev
RUN apt-get install -y libsqlite3-dev
RUN apt-get install -y gdb
RUN apt-get install -y git

RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 1

# set environment variables
ENV HOME=/home

# --------- DPMC ---------
ENV DMPC=$HOME/dpmc
ENV ADDMC=$DPMC/addmc
ENV DMC=$DPMC/dmc
ENV ADDMC_COLAMD=$ADDMC/libraries/colamd
ENV ADDMC_LIBS=$ADDMC/libraries
ENV ADDMC_SYLVAN=$ADDMC/libraries/sylvan/
ENV SUITE_SPARSE=$ADDMC_LIBS/SuiteSparse
ENV DMPC_WEIGHTED_TESTS=$DPMC/tests/weighted
ENV DMPC_UNWEIGHTED_TESTS=$DPMC/tests/unweighted
ENV LG=$DPMC/lg
ENV LG_SOLVERS=$LG/solvers
ENV HTD_SOLVER=$LG_SOLVERS/htd-master
ENV DPMC_BUILD=$HOME/build/dpmc
RUN cd $BUILD && mkdir dpmc


# compile the colamd library and copy it into the addmc library directory
RUN cd $SUITE_SPARSE/SuiteSparse_config/ && make -s
RUN rm -f $ADDMC_COLAMD/libcolamd.a
RUN cd $SUITE_SPARSE/COLAMD/ && make -s
RUN cp $SUITE_SPARSE/COLAMD/Lib/libcolamd.a $ADDMC_COLAMD/

# prepare the dmc executable
RUN cd $DMC/ && make dmc -j 8
RUN cp $DMC/dmc $DPMC_BUILD

# prepare the lg executable
RUN cd $LG/ && make -j 8
RUN cp $LG/build/lg $DPMC_BUILD

# prepare the htd_main-1.2.0 executable
RUN cd $HTD_SOLVER && mkdir build
RUN cd $HTD_SOLVER/build && cmake .. && make -j 8
RUN cp $HTD_SOLVER/build/bin/htd_main-1.2.0 $DPMC_BUILD

RUN cp $DMPC_WEIGHTED_TESTS/mcc21__wff.3.75.315.cnf $DPMC_BUILD


