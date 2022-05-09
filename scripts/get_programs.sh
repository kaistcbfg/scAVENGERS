#!/bin/bash

wget -O $1/strelka-2.9.10.centos6_x86_64.tar.bz2 -P $1 \
https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
wget -O $1/2.0.tar.gz -P $1 \
https://github.com/wheaton5/souporcell/archive/refs/tags/2.0.tar.gz
tar -xvf $1/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xvzf $1/2.0.tar.gz
cd $1/souporcell-2.0/troublet && cargo build --release
mv $1/souporcell-2.0/troublet/target/release/troublet $1/troublet
ls $1/strelka-2.9.10.centos6_x86_64/bin/
rm -r $1/souporcell-2.0 $1/2.0.tar.gz $1/strelka-2.9.10.centos6_x86_64.tar.bz2