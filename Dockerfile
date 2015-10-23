FROM andrewosh/binder-base

MAINTAINER Zach Sailer <zachsailer@gmail.com>

USER root

# Download and build cdhit
RUN cd $HOME

RUN git clone git@github.com:weizhongli/cdhit.git

RUN cd cdhit
RUN ./configure
RUN make
RUN mv cd-hit cdhit

RUN cd $HOME

# download and build phyml
RUN git clone git@github.com:stephaneguindon/phyml.git

RUN cd phyml
RUN ./configure
RUN make

RUN cd $HOME

# Download and build
RUN wget http://sourceforge.net/projects/msaprobs/files/MSAProbs-0.9.7.tar.gz > msaprobs.tar.gz
RUN tar -xzvf msaprobs.tar.gz

RUN cd $HOME

# Export the path to these commands
RUN PATH=$HOME:$PATH/bin
RUN export PATH
RUN export PATH=$HOME/cdhit:$PATH
RUN export PATH=$HOME/msaprobs:$PATH
RUN export PATH=$HOME/phyml-20120412:$PATH

RUN pip install phylogenetics
