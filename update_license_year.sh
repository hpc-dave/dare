#!/bin/bash

find -type f -not -path '*/\.*' -exec sed -i 's/(c) 2023/(c) 2023/g' {} +
