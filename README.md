# iso-filt
Filtration of RNA isoforms for selected genes from long-read RNA-seq data.

Alignment of long-reads from RNA-seq often results in various levels of noise due to  factors such as DNA-contamination, PCR-artifacts and misalignments, which may challenge downstream interpretations of an experiment. Furthermore, downstream analyses may have limited capacity to handle large BAM files.

To handle these challenges, this script filters reads for selected genes so that only those reads that match filtering criteria get included.

The included reads are written to a new, sorted BAM file that only contains the filtered reads for the selected genes.

Two filtration methods are used and each generates a BAM file:

- Full length filtering
- Isoform filtering

### Full length filtering
Includes reads that span the full transcript length and overlap the first and the last exons. This method may not be very useful for long genes and it may be less efficient at filtering out genomic contamination.

### Isoform filtering
Includes reads that match at least two annotated exons. This method typically generates the most reads because it includes both full length reads and reads that do not
span the full transcript length. The latter is quite frequent with methods based on cDNA and PCR because of RNA breakdown and internal priming.

In addition, this method tends to generate cleaner results because an "exon match" is more restrictive than an overlap.

### Requirements
- Pysam
- Annotation file matching the reference used for alignment

### Running iso-filt
For a small help text explaining how to run iso-filt, run the following command in a directory containing iso_filt.py:
´´´
python iso_filt.py -h
´´´

### Important
This script is designed for long-read data and will most likely not work with short-read data.
