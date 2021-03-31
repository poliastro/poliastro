#!/usr/bin/bash

pep8 .
pytest --cov-report term-missing --cov . src/
pylint poliastro --disable=no-member,no-name-in-module,invalid-name --reports=n
