FROM ubuntu:22.04

WORKDIR /docker_build/

# Install required packages
RUN apt-get update
RUN apt-get install -y build-essential libbz2-dev libcurl4-openssl-dev libssl-dev zlib1g-dev liblzma-dev libdeflate-dev libncurses5-dev curl git

RUN git clone https://github.com/e-jorsboe/fastNGSadmix.git
RUN cd fastNGSadmix && make && cp fastNGSadmix /usr/local/bin

WORKDIR /
