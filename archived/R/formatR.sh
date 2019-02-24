#!/bin/bash

#convert to format that both unix and windows can read
find ./ -iname '*.R' -exec sed -i 's,\.\\\\,,g' '{}' ';'
find ./ -iname '*.R' -exec sed -i 's,\\\\,\/,g' '{}' ';'
find ./ -iname '*.R' -exec sed -i 's,\/\/,\/,g' '{}' ';'

# replace the <- which is just ugly as so many things in R
find ./ -iname '*.R' -exec sed -i 's/<-/=/g' '{}' ';'
