#!/usr/bin/env bash

SUPPORTED_VERSIONS=(33 34 35)
BUILD_NUMBER=1

python setup.py sdist
python setup.py register
python setup.py sdist upload

for ii in $SUPPORTED_VERSIONS
    do conda build . --python $ii
done
for ii in $SUPPORTED_VERSIONS
do conda convert \
    /home/juanlu/.miniconda3/conda-bld/linux-64/poliastro-0.4.0-py${ii}_${BUILD_NUMBER}.tar.bz2 \
    --platform all
done
find . -name "*.tar.bz2" | xargs -P0 -I{} binstar upload {} -u poliastro
