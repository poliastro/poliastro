DOCKER_BASE_IMAGE := "poliastro:dev"
DOCKER_CONTAINER_NAME := "poliastro-dev"


image: Dockerfile pyproject.toml
	docker build \
	  -t poliastro:dev \
	  .

docker:
	docker run \
	  -it \
	  --rm \
	  --name ${DOCKER_CONTAINER_NAME} \
	  --volume $(shell pwd):/code \
	  --user $(shell id -u):$(shell id -g) \
	  ${DOCKER_BASE_IMAGE} \
          bash

.PHONY: docker image
