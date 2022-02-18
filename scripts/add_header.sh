#!/bin/bash

FIRST=$(head -1 $1)
echo $FIRST
if [[ $FIRST != "x"* ]]; then
    echo "no header"
    sed -i '1s/^/x\ty\tz\tgene\n/' $1
fi
