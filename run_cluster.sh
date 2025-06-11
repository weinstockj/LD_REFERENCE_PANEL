#!/bin/bash

snakemake -p --profile slurm -j 500 1>stdout.log 2>stderr.log
