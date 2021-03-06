#!/bin/bash

file="include/version.hpp"

echo "// File generated by writeVersion.bat/sh. Do not modify" > $file

echo -n "#define GIT_HEAD \"" >> $file
git log --pretty=format:'%h' -n 1 >> $file
echo "\"" >> $file
echo "" >> $file

echo -n "#define GIT_LOG " >> $file
logOutput=$(git log -1)
logOutputNew=$(echo "$logOutput" | sed -e 's/\\/\\\\/g;s/"/\\\"/g;s/((\\)+"*)*/\\1/g;s/\r//g;s|$|\\n" \\|g;s/^/"/g')

echo "$logOutputNew" >> $file

echo "\"\n\"" >> $file
echo "" >> $file

echo -n "#define GIT_DIFF " >> $file
diffOutput=$(git diff)
diffOutputNew=$(echo "$diffOutput" | sed -e 's/\\/\\\\/g;s/"/\\\"/g;s/((\\)+"*)*/\\1/g;s/\r//g;s|$|\\n" \\|g;s/^/"/g')

echo "$diffOutputNew" >> $file

echo "\"\n\"" >> $file
echo "" >> $file

exit 0