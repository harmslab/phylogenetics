FROM andrewosh/binder-base

MAINTAINER Zach Sailer <zachsailer@gmail.com>

USER root

# Download and build cdhit
RUN cd $HOME

RUN git clone https://github.com/weizhongli/cdhit.git

RUN cd cdhit \
    ./configure \
    make \
    mv cd-hit cdhit

# download and build phyml
RUN git clone https://github.com/stephaneguindon/phyml.git

RUN  cd phyml \
    ./configure \
    make

RUN cd $HOME

# Download and build
RUN curl -L http://sourceforge.net/projects/msaprobs/files/latest/download > msaprobs.tar.gz
RUN tar -xzvf msaprobs.tar.gz

RUN cd $HOME

# Export the path to these commands
RUN PATH=$HOME:$PATH/bin
RUN export PATH
RUN export PATH=$HOME/cdhit:$PATH
RUN export PATH=$HOME/msaprobs:$PATH
RUN export PATH=$HOME/phyml-20120412:$PATH

# clone and install phylogenetics
RUN git clone https://github.com/Zsailer/phylogenetics

RUN cd phylogenetics \
    python3 setup.py install \
    python2 setup.py install
