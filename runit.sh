#!/bin/sh
if [ ! -d "build" ]; then
  ./buildit.sh
fi
build/bin/$1 scenes/$2.json
