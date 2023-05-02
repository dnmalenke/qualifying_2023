#! /bin/bash

if [ -d _build ]; then
    cd _build
else
    mkdir _build
    cd _build
    cmake -DCMAKE_BUILD_TYPE=Debug ..
fi

if ! cmake --build .; then
    exit 1;
fi

cd bin
if  ./EC2023; then
    cat ./scores.json
fi
