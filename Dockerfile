FROM python:3.12.9

COPY . /ipensive
WORKDIR /ipensive
# RUN rm /ipensive/config.py

RUN git config --global http.sslVerify false
RUN pip install -e .