LABEL "repository"="https://github.com/fcakyon/conda-publish-action"
LABEL "maintainer"="Fatih C Akyon"

FROM ubuntu:16.04
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt-get update
RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

# install miniconda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
RUN conda --version


RUN apt-get update
# to fix: import cv2 > ImportError: libGL.so.1: cannot open shared object file: No such file or directory
RUN apt-get install -y libgl1-mesa-dev

# to fix: import cv2 > ImportError: libjasper.so.1: cannot open shared object file: No such file or directory
RUN apt-get install -y libjasper1

RUN conda install -y anaconda-client conda-build conda-verify

COPY entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]