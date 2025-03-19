FROM python:3.7.13

COPY . /ipensive
WORKDIR /ipensive

RUN pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org -r requirements.txt
