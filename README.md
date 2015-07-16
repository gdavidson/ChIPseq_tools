#Stand-alone scripts to work on ChIP-seq data.

## getFromBed.py

```
Allows operations on BED Files.
Arguments:
 -i <file1.bed>, input file in bed format
Options (choose one):
 -o <offset(integer)>, changes the coordinates of each entry of the input file by the offset amount (start-offset, end+offset)
 -b <file2.bed>, finds all entries in file2 which ID/name matches any ID/name in file1 (4th column) .
Optional:
 -f <field number(integer)>, only with '-b', field (column) which contains IDs for file2.bed (default:4)
```

## getFromFasta.py

```
 Allows different operations on fasta files.
 Arguments:
 -i <file.fasta>, sequence(s) file in fasta format
 Options (choose one):
 -r <regexp>,  looks for the regular expression within input sequences.
 -u true, generates a new multi fasta file with new uniques IDs for each sequence.
 -f <ids.txt>, file containing a list of sequences IDs (1 per line) to retrieve from '-i' file.
 Optional:
 -p <separator>, only if '-f', parses the IDs of both files (ex: sequenceID(file.fasta): MITF_peak_201, ID(ids.txt): CT_HG19_MITF_peak_201, then use either 'MITF' or 'peak' as separator).
```

## getFromAnnotations.py

```
 Allows operations on homer 'annotatePeaks.pl' output file.
 Arguments:
 -i <homer_annotation.txt>, homer output file.
 Options:
 -h <histogram title>, makes a histogram of distances to nearest TSS.
 -p <piechart title>, makes a chart summarizing annotations (number of peaks annotated as 'promoter-TSS', 'exon, 'intron' ect...).
 -r <keyword>, retrieves lines annotated as the keyword (<promoter-tss>, <exon>, <intron>, <TTS>, <intergenic>, <non-coding>, <3'UTR>, <5'UTR>)
 -b <biomart_export.txt>, completes annotation with a biomart output file containing the following fields (in that order): transcript ID, gene ID, gene name, description. Adds these three fields to the homer file.
```

## editGTFlocations.py

```
Changes first field of a GTF file to match UCSC standards (chr1, chr2 ...).
Arguments:
 -i <file.gtf>, the GTF file to edit.
 -a <assembly>, name of the assembly (hg19 or mm9).
Optional:
 -p <mm9_patches.tsv>, the file containing patches and loci for your assembly (alt_scaffold_placement.txt in genebank genome DB).
```

## compareGeneLists.py

```
 Compares a gene list to a homer annotatePeaks.py output and finds common genes.
 Arguments:
-i <genes.txt>, gene list with two tab separated fields per line (ENSEMBL Gene ID, gene common name).
-h <homerOut.tsv>, annotated file from homer with ensembl transcriptID as reference.
-e <biomartOut.tsv>, biomart output for your ensembl release and organism (two columns: geneID, transcriptID). If -g, first column has to be geneID then other fields are optional. 
 Optional:
-b true, if your '-i' file is also a homer annotation file.
-g true, if your '-h' file is also a list with gene IDs in the first column.  
```

## clinvarToBed_iterative.py

```
 Converts the ClinVar DB XML file to bed format
 Usage: python clinvarToBed_iterative.py -i ClinVarFullRelease_2015-04.xml
 Args:
-i: input file, clinVar XML file
 Output:
'clinvar.bed': clinVar DB in bed format
```

## getMotifLocations.py

```
Writes the genomic locations of a MEME-found motif in bed format.
 Arguments:
-m <meme_motif.txt>, copy paste the locations from the MEME HTML output into a .txt file (example line:'1684.    hg19_ct_macs2summits_2357_sox10_peak_3016    -    49    1.24e-8    ACAACAACAC    AAAAGGCCCCTTTGT    TACGGCCCTG')
-p <macs_peaks.bed>, MACS predicted peaks in bed format.
-w <motif width(integer)>, width of the motif (you can fiddle with the size)
-n <experiment name>, name of the experiment, must be the '-n' argument you used in MACS.
 Optional:
-f true, if your '-m' file is an output from FiMo
```

## getSNVTable.py

```
 Writes a table summarizing single-nucleotide variants (SNVs) found in a motif.
 Arguments:
 -i <intersect_motif_peaks.bed>, output file from intersectBed with '-a': clinvar Bed file, '-b': the motif locations (from getMotifLocations.py) and '-wo'.
 Optional:
 -v <clinvarMain.bed>, clinvar bed file (with all columns) containing more information on the variants. The output will contain more details on each SNV.
 -h <homer_annotation.txt>, homer annotation  of the peak file containing the motif. The output will contain annotation details of the peak containing the variant.
 -d <disease_name.txt>, clinvar disease file (ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/disease_names)
```

