# TraceTrack

Sanger sequencing analysis web app. View the [publication](https://biorxiv.org/cgi/content/short/2022.07.28.501824v1) or try the app [online](https://tracetrack.dichlab.org).

## Overview
The TraceTrack application serves for automatic batch analysis of Sanger sequencing trace files. It
accepts a set of trace files and a reference file, performs alignments of trace sequences to 
corresponding reference sequences and displays an overview of resulting alignments.

For more information on how to use the app, see the "Help" page in the application itself.


## Run TraceTrack in Docker container

### 1. Install Docker

See https://docs.docker.com/docker-for-mac/install/

### 2. Build all images using Docker Compose

```bash
make docker-build
```

### 3. Run all services using Docker Compose

```bash
make docker-run
```

To build and run, you can use:
```bash
make docker-build docker-run
```
### 4. Handle code updates

After your code is updated, you will need to stop the services, run build and start again. 
See the next section for info on running locally with flask auto-reload.

### 5. Frontend

Frontend is available through your browser at http://localhost:5001/

## Run TraceTrack using Conda

### 1. Install Conda

See https://docs.conda.io/en/latest/miniconda.html

### 2. Install Redis server

You can [install Redis using Brew](https://medium.com/@petehouston/install-and-config-redis-on-mac-os-x-via-homebrew-eb8df9a4f298).

### 3. Setup environment

```bash
# Install dependencies using environment.yml
make install
# Update environment with:
make update-env
```

### 4. Run all services

You will have to run each service in a separate terminal (Use Cmd+T to open a new tab):

```bash
# Run Redis server
make redis

# In a separate terminal, run celery worker queue
make celery

# In a separate terminal, run flask web server
make web
```
### 5. Handle code updates

After your code is updated, the flask web service should refresh automatically. However, the celery service
needs to be stopped and started manually, so you will need to do that if you update code 
that is executed from the workers.

### 6. Frontend

Frontend is available through your browser at http://localhost:5001/

## Files

- [Dockerfile](Dockerfile) 
- [Makefile](Makefile) 
- [README.md](README.md) 
- [CHANGELOG.md](CHANGELOG.md) Descriptions of release versions
- [tracetrack](tracetrack) 
  - [flask_server.py](tracetrack/flask_server.py) Basic app functionality, processing uploaded files
  - [server_utils.py](tracetrack/server_utils.py) Functionality for displaying, downloading and exporting data
  - [alignment_utils.py](tracetrack/alignment_utils.py) Contains functions for sequence processing and alignment using Clustal Omega
  - [tasks.py](tracetrack/tasks.py) Task for task queue
  - [celeryconfig.py](tracetrack/celeryconfig.py) Configuration of task queue
  - [loading_data.py](tracetrack/loading_data.py) Loading and processing input files functionality
  - [web.py](tracetrack/web.py) CLI for running with no task queue
  - [entities](tracetrack/entities)
    - [alignment.py](tracetrack/entities/alignment.py) Contains Alignment class
    - [aligned.py](tracetrack/entities/aligned.py) Contains AlignedTrace and AlignedReference classes
    - [errors.py](tracetrack/entities/errors.py) Contains exception classes
    - [record.py](tracetrack/entities/record.py) Contains classes Feature, Record, TraceSeqRecord and find_quality_ends function
    - [reference_db.py](tracetrack/entities/reference_db.py) Contains ReferenceDb class and methods from reading reference sequences from files
    - [sequence_position.py](tracetrack/entities/sequence_position.py) Contains SequencePosition class
  - [templates](tracetrack/templates) 
    - [index.html](tracetrack/templates/index.html)  Source for main page of application, contains functionality for uploading files 
    - [layout.html](tracetrack/templates/layout.html) General layout of page
    - [help.html](tracetrack/templates/help.html) Help page
    - [loading.html](tracetrack/templates/loading.html) Page displayed while loading results
    - [results.html](tracetrack/templates/results.html) Source for page for displaying results table with alignments
    - [settings.html](tracetrack/templates/settings.html) Source of page for choosing alignment settings, such as quality trimming
    - [welcome.htlm](tracetrack/templates/welcome.html) Source for welcome page
  - [static](static)
    - [stylesheet.css](stylesheet.css) Defines style classes for use in the templates
    - [trace_viewer.js](tracetrack/static/trace_viewer.js) Script for showing trace file chromatograms
    - [table_display.js](tracetrack/static/table_display.js) Javascript functions for setting up results table
- [data](data)
   - [example1](example1) Contains example trace files and reference sequence tables
- [tests](tests) 
  - [test_alignment.py](tests/test_alignment.py) Tests for alignment module
  - [test_loading_data.py](tests/test_loading_data.py) Tests for loading_data module
- [resources](resources) Data for tests
- [environment.yml](environment.yml)
- [docker-compose.yml](docker-compose.yml) Configuration of task queue container
- [setup.py](setup.py) Script for tracetrack package installation
- [jenkins](jenkins)
