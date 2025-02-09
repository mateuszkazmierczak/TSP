#!/bin/bash


files=(
  bier127
  
  # ... 
)

for f in "${files[@]}"; do
    echo "Running ./main samples/$f.txt parameters/$f.txt ..."
    ./main samples/"$f".txt parameters/"P_berlin52".txt
done
