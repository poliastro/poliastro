#!/usr/bin/bash

pep8 .
py.test-3.3 --cov-report term-missing --cov . tests/
pylint poliastro --disable=no-member,no-name-in-module,invalid-name --reports=n
