#!/bin/bash

snakemake -p --profile slurm -j 50 1>stdout.log 2>stderr.log
