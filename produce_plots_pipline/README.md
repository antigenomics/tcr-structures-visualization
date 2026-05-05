# Produce Plots Pipeline

This folder contains the pipeline code for producing plots and processing structure data.

- `alignment.py` - alignment utilities
- `coordinates.py` - coordinate handling utilities
- `plotting.py` - plotting helpers
- `run_pipeline.py` - pipeline orchestration
- `get_files.sh` - helper script for file management

To run pipline use python3 -i /folder/with/pdb/files. Strucutres should be described in /projects/structures/clusters/HomoSapiens_MHCI_all_clusters.tsv or simiilar files saved to `data_generation` variable.
