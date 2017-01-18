#!/bin/sh
# install scripts and configuration needed for HSRL archive computers
prefix=/usr/local/hsrl
cp -pr lib ${prefix}
cp -pr config ${prefix}
cp -pr etc ${prefix}/
cp -pr bin /usr/local/

mkdir -p ~eol-lidar/etc
cp -p eol-lidar/etc/master.crontab ~eol-lidar/etc
