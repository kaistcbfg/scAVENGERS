#!/bin/bash

wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
wget https://github.com/wheaton5/souporcell/archive/refs/tags/2.0.tar.gz
tar -xvzf strelka-2.9.10.centos6_x86_64.tar.gz
tar -xvzf 2.0.tar.gz
mv souporcell-2.0/troublet .
cd troublet && cargo build --release
rm -r souporcell-2.0 2.0.tar.gz strelka-2.9.10.centos6_x86_64.tar.gz