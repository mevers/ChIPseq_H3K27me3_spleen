#!/bin/bash

src=$1
dir=$2

rsync -avzh --progress $src $dir
