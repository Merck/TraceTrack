FROM continuumio/miniconda3

RUN apt-get update -y \
        && apt-get install -y --no-install-recommends \
            libgtextutils-dev libtbb-dev \
            autoconf automake libcurl4-gnutls-dev libncurses5-dev \
            build-essential pkg-config \
        && /usr/local/tracetrack/clean-apt-and-logs.sh

ENV PATH /opt/miniconda/tracetrack:$PATH
WORKDIR /tmp

COPY environment.yml environment.yml
RUN conda env update -n base -f environment.yml \
    && rm -rf /tmp/* \
    && conda clean --all \
    && rm -r /opt/miniconda/lib/python3.7/site-packages/future/backports/test/ /opt/miniconda/pkgs/*

WORKDIR /opt/suprseqr/

COPY tracetrack tracetrack
COPY data data
COPY tests tests
COPY resources resources
COPY setup.py setup.py

RUN pip install .

RUN useradd docker \
  && mkdir /home/docker \
  && chown docker:docker /home/docker \
  && addgroup docker staff
USER docker

ENV FLASK_APP tracetrack.flask_server
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD [ "flask", "run", "--host", "0.0.0.0" ]





