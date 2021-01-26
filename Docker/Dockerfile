FROM nlknguyen/alpine-mpich
RUN sudo apk add --no-cache zlib-dev git
RUN git clone https://github.com/LANL-Bioinformatics/PCRamp.git pcramp \
    && cd pcramp \ 
    && make \
    && sudo cp pcramp /usr/local/bin/ \
    && cd .. && rm -rf pcramp

CMD ["pcramp"]

## To build: docker build -t pcramp .
## Testing/development Usage:  docker run -it --rm -v $(pwd):/project
## Docker compose (Single Host Orchestration), Docker Swarm (Multi-Host Orchestration)
##            Usage see: https://github.com/NLKNguyen/alpine-mpich/tree/master/cluster