# Docker image for poliastro development
FROM python:3.7
LABEL maintainer="Juan Luis Cano Rodríguez <hello@juanlu.space>"

# https://pythonspeed.com/articles/activate-virtualenv-dockerfile/
ENV VIRTUAL_ENV=/opt/venv
RUN python -m venv ${VIRTUAL_ENV}
ENV PATH="${VIRTUAL_ENV}/bin:${PATH}"

RUN python -m pip install -U pip setuptools
RUN python -m pip install flit pygments wheel

WORKDIR /code

# Unfortunately this will invalidate the cache
# for every code change
# Wouldn't it be nice to have a flit compile
# that freezes a requirements.txt?
COPY . /code

ENV FLIT_ROOT_INSTALL=1
RUN flit install --symlink --extras dev
