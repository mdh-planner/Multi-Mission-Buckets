# Multi-Robot Allocation and Optimization in a Multi-Mission Framework

This repository contains the experimental code and scripts used for our IFAC paper on **multi-mission multi-robot planning with time-bucket reservations and replanning**.

The C++ code builds and solves mixed-integer optimization models (CPLEX) for:
- A **first plan** over a finite horizon.
- A **cut & replan** phase after a chosen cut time, with **pseudo-depots**, **robot release times**, and optional **new tasks**.
- Multiple baselines: **OURS (reservations)**, **no bucket layer**, **monolithic**, etc.

MATLAB scripts post-process the exported CSV/JSON files to reproduce **Gantt charts** and **bucket-reservation diagrams** like those in the paper.
