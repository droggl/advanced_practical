#!/bin/bash

# Compile your program, ie:
make

# Check compilation succeeded
if [[ $? -ne 0 ]]; then
    echo "compilation failed"
    exit 1
fi

# Pipe your input into the application
./application $1 $2
