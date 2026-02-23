#!/bin/bash
# Simplified EGAPx runner - container has singularity installed internally

# Make cache directory if it doesn't exist
mkdir --parents /tmp/${USER}/singularity_cache

# Set NXF_SINGULARITY_CACHEDIR to the cache directory
export NXF_SINGULARITY_CACHEDIR=/tmp/${USER}/singularity_cache

# Run the EGAPx container with the provided arguments
singularity run --cleanenv --bind /mmfs1/:/mmfs1/ srlab-NCBI-EGAPx.sif "$@"
