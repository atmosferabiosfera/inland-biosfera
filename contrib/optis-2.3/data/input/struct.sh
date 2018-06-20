#!/bin/bash

mkdir -p data/offline/single_point/;

mkdir data/offline/single_point/output;
mkdir data/offline/single_point/params;

ln -s data/offline/single_point/conf conf;
ln -s data/offline/single_point/params params;
ln -s data/offline/single_point/input input;
ln -s data/offline/single_point/output output;

ln -s /home/user/Programas/inland/inland-single_point inland;

cp /home/user/Programas/optis/data/input/params/* params;