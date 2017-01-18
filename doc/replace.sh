#!/bin/bash

for f in source/${1}*.rst ; do
		mv $f $( echo $f | sed s,/${1}.,/,g )
done