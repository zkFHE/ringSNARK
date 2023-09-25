FROM ubuntu:22.04
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y cmake clang git libboost-all-dev
COPY . /src

WORKDIR /src
RUN mkdir -p build_docker
WORKDIR /src/build_docker
RUN cmake ..
RUN cmake --build .

CMD ./bin/bench_logistic_regression_inference