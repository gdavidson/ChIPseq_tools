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






