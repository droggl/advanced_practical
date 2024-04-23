#!/bin/bash

# Compile your program, ie:
clang++ -std=c++17  $1 -o application

# Check compilation succeeded
if [[ $? -ne 0 ]]; then
    echo "compilation failed"
    exit 1
fi

# Pipe your input into the application
./application $2 $3
