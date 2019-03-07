#!/bin/bash
cat $1 | perl -pe 'BEGIN { binmode STDIN, ":utf8"; } s/(.)/ord($1) < 128 ? $1 : sprintf("\\U%08x", ord($1))/ge;' > /tmp/$1
