#!/bin/sh

# Bash script for automatically retrieving all the Roberts Lab repos
# and their corresponding wikis.

mkdir ~/gitrepos

mkdir ~/gitrepos/RobertsLab

cd ~/gitrepos/RobertsLab

git clone git@github.com:RobertsLab/onboarding.git


git clone git@github.com:RobertsLab/onboarding.wiki.git

git clone git@github.com:RobertsLab/resources.git


git clone git@github.com:RobertsLab/resources.wiki.git


git clone git@github.com:RobertsLab/code.git

git clone git@github.com:RobertsLab/code.wiki.git
