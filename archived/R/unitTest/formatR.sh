#!/bin/bash

#convert to format that both unix and windows can read
sed -i 's,\.\\\\,,g' *.R
sed -i 's,\\\\,\/,g' *.R
sed -i 's,\/\/,\/,g' *.R

# replace the <- which is just ugly as so many things in R
sed -i 's/<-/=/g' *.R
