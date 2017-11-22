#same as the build host !
#FROM centos:7.2.1511 

#alternate slim images
#FROM bitnami/minideb:jessie
FROM frolvlad/alpine-glibc:alpine-3.5

WORKDIR /app 
COPY bin/ bin/
COPY lib/ lib/

#CMD ["bin/infoli.x.docker","-n","100"]
