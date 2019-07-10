FROM ubuntu:latest

WORKDIR /build

COPY . .

RUN apt-get update -qq && \
    apt-get install -y zlib1g-dev \
                    git \
                    cmake \
                    build-essential \
                    python3
RUN cd Big-BWT; make; cd ..; 
ENV PATH /build/Big-BWT:$PATH
RUN rm -rf build; mkdir -p build; cd build; cmake ..; make; make install

WORKDIR /app

CMD ['bash']
