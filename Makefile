IMAGE_NAME := tracetrack
IMAGE_URL := $(IMAGE_NAME)_web

docker-build:
	docker-compose build

docker-terminal:
	docker run -v --rm -it $(IMAGE_URL) /bin/bash

docker-run:
	docker-compose up

docker-test:
	docker run --rm $(IMAGE_URL) pytest tests

#-------------------#
# Local environment #
#-------------------#

CONDA_ACTIVATE=source $$(conda info --base)/bin/activate
ACTIVATE_LOCAL_ENV= $(CONDA_ACTIVATE) tracetrack-local; \
					export CELERY_BROKER_URL="redis://localhost:6379/0"; \
					export CELERY_RESULT_BACKEND="redis://localhost:6379/0"

update-env: environment.yml
	conda env update -n tracetrack-local -f $< && \
	$(ACTIVATE_LOCAL_ENV) && \
	pip install -e .

install:
	conda create -n tracetrack-local && \
	$(ACTIVATE_LOCAL_ENV);
	@$(MAKE) update-env

redis:
	$(ACTIVATE_LOCAL_ENV); \
	redis-server

celery:
	$(ACTIVATE_LOCAL_ENV); \
	C_FORCE_ROOT=1 celery -A tracetrack.tasks worker --loglevel=info

flower:
	$(ACTIVATE_LOCAL_ENV); \
	flower -A tracetrack.tasks --port=5555

web:
	$(ACTIVATE_LOCAL_ENV); \
	export FLASK_APP=tracetrack.flask_server; \
	export FLASK_DEBUG=1; \
	flask run --host 0.0.0.0 -p 5001
