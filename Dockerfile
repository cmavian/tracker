#docker pull the latest image of nextflow, nextclade
FROM nextstrain/nextclade:latest as nextclade
FROM nextflow/nextflow:22.10.7 as nextflow
# merge
COPY --from=nextclade ./usr/bin/nextclade /usr/bin/nextclade
RUN /bin/sh -c set -eux && ln -s /usr/bin/nextclade /nextclade

# install python3.9
RUN yum install -y epel-release gcc openssl-devel bzip2-devel libffi-devel zlib-devel tar vim nano
RUN yum groupinstall -y "Development Tools"
RUN yum install -y python-pip
RUN curl -L -o 'Python-3.9.6.tgz' 'https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz'
RUN tar -xvf Python-3.9.6.tgz
RUN rm Python-3.9.6.tgz
WORKDIR /Python-3.9.6
RUN ./configure --enable-optimizations
RUN make altinstall
ENV PYTHON_VERSION=3.9.6
RUN /usr/local/bin/python3.9 -m pip install --upgrade pip
RUN pip install biopython

#specify the manteiner of the docker
MAINTAINER Carla Mavian <cmavian@ufl.com>  Devon Gregory <virologistd@gmail.com>

#setting the current work directory
WORKDIR /work
RUN mkdir data


#python scripts input
ADD ./bin /work/bin
ADD data/mutations_list.txt /work/data
ADD data/*.fa /work/data

CMD ["bash"]

# execute first time workflow with : bash bin/tracker.sh -m /work/data/mutations_list.txt data/test.fa
# execute then workflow with : bash bin/tracker.sh data/test.fa
