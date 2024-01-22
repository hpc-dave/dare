#!/bin/bash

find -type f -not -path '*/\.*' -exec sed -i 's/(c) 2024/(c) 2024/g' {} +
