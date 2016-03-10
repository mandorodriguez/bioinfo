## Changes to xml 

We removed the sam to bam conversion after the Bowtie2 run, due to memory issues with large batches of data. The sam/bam conversion is
done in a second step in the workflow. We also modified the naming of the output file to propagate the name from the input files to the 
output. Galaxy workflow renaming of output files does not work with lists as input. 
