### Usage

1. Run `load_data.R`:

`Rscript load_data.R data doublets out`

where `data` is a folder with samples (each sample in its own folder), `doublets` is a file with a single column where 0-s and 1-s correspond to whether a barcode was identified as a singlet (0) or a multiplet (1), `out` - output folder.

2. Run `sclean.R`:

`Rscript sclean.R data out`

where `data` is the pre-processed object-list (output file of load_data.R), `out` - output folder.

3. Run `integrate.R`:

`Rscript integrate.R data out`

where `data` is the pre-processed and QC'ed object-list (output file of sclean.R), `out` - output folder.