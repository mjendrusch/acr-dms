#!/bin/bash


if conda info --envs | grep acrdms > /dev/null; then
  echo "already set up."
else
  echo "setting up..."
  # create conda environment
  conda create -n acrdms -y python=3.6

  # activate conda env
  conda activate acrdms

  # install packages
  conda install pytorch cpuonly -c pytorch -y
  conda install -c conda-forge -y notebook matplotlib numpy scikit-learn biopython seaborn
  pip install -q FlowCal
  conda deactivate
  echo "done!"
fi

conda activate acrdms
