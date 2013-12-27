#!/usr/bin/bash

pep8 .
py.test --cov-report term-missing --cov . tests/
pylint poliastro --disable=no-member,no-name-in-module,invalid-name
