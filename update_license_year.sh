#!/bin/bash

find -type f -not -path '*/\.*' -exec sed -i 's/(c) 2025/(c) 2025/g' {} +
