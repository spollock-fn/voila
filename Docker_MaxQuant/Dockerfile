FROM ubuntu:20.04

RUN apt-get -y update
RUN apt-get -y install wget vim unzip

RUN wget https://packages.microsoft.com/config/ubuntu/20.04/packages-microsoft-prod.deb -O /tmp/packages-microsoft-prod.deb
RUN dpkg -i /tmp/packages-microsoft-prod.deb
RUN rm /tmp/packages-microsoft-prod.deb

RUN apt-get update && \
    apt-get install -y apt-transport-https && \
    apt-get update && \
    apt-get install -y aspnetcore-runtime-2.1

RUN apt-get update && \
    apt-get install -y apt-transport-https && \
    apt-get update && \
    apt-get install -y dotnet-sdk-2.1

RUN apt-get -y clean

COPY MaxQuant /MaxQuant
#RUN unzip max_quant.zip && rm max_quant.zip

WORKDIR /MaxQuant

# COPY mqpar.xml test.raw uniprot.fasta ./
# mono MaxQuant/bin/MaxQuantCmd.exe mqpar.xml -n

CMD [ "/usr/bin/dotnet", "bin/MaxQuantCmd.exe", "data/mqpar.xml" ]
#CMD [ "/usr/bin/dotnet", "--version" ]
