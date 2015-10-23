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

# Download sip
#RUN curl -L http://www.riverbankcomputing.com/software/sip/download > sip.tar.gz
#RUN tar -xzvf sip.tar.gz
#RUN cd sip \
#    python configure.py \
#    make \
#    make install \

# Download pyqt4
#RUN curl -L http://www.riverbankcomputing.com/software/pyqt/download >pyqt.tar.gz
#RUN tar -xzvf pyqt.tar.gz
#RUN cd pyqt \
#    python configure.py \
#    make \
#    make install \

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

# Install phylogenetics
RUN pip install phylogenetics ete2 numpy
