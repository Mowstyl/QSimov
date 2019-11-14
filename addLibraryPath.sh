#!/bin/bash

# We add the lib dir to the ldd search path
echo "$(pwd)/lib" > /etc/ld.so.conf.d/libqsimov.conf
# We force ldd to reload the conf files
ldconfig
