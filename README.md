# splitStrains

### TODO:
1) Need to run a check on the gff file (version).  
2) Separate *splitStrains* code into files.

### Notes:    
Make sure that <b>sorted and indexed BAM's</b>  aligned sequences and <b>indexed </b> fasta reference have the same sequence ID.  
For example, if BAM aligned sequences refer to "*gi|41353971|emb|AL123456.2|* " then the fasta reference file should start with "*>gi|41353971|emb|AL123456.2|* ".

### Possible bad behavior:
When depth coverage is high (300 and greater) it is possibe that likelyhood_ratio_test function can overflow.  
Couldn't reproduce the error so far.  
