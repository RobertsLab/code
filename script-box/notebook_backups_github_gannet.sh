#!/bin/bash

# Shell script to backup lab notebooks hosted on GitHub.

# Notebook path
notebook_dir="/volume2/web/github_backups/notebooks"

# Array of notebooks
notebooks_array=("${notebook_dir}/grace-ac.github.io" "${notebook_dir}/LabNotebook" "${notebook_dir}/sams-notebook" "${notebook_dir}/shellytrigg.github.io" "${notebook_dir}/sr320.github.io" "${notebook_dir}/yaaminiv.github.io")

# Loop through notebook dirs and pull any repo updates
for notebook in "${notebooks_array[@]}"
do
	cd "${notebook}" || exit
	git pull
	cd "${notebook_dir}" || exit
done
