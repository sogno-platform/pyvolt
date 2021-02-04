FROM python:3.7

ADD . /pyvolt
WORKDIR /pyvolt

RUN python3 setup.py develop