#!/usr/bin/bash

pep8 .
py.test --cov . tests/
pylint poliastro --disable=no-member,no-name-in-module,invalid-name
