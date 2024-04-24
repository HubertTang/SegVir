## Missing segment datasets

`missing_seg.SRR16905247.fna` and `missing_seg.SRR19790906.fna` contain the non-RdRp segmented RNA viral contigs identified from the missing segment datasets. Considering the possibility of false positives in SegVir’s output, we aligned the contigs identified by SegVir with the NCBI’s nt database and removed several contigs that may originate from eukaryotic or prokaryotic organisms. The contigs shown in the files are after this filtering process. In all samples, only a single virus’s RdRp was found under the same family in each sample. Therefore, we believe that the non-RdRp contigs with the same function found under this family all come from the same segment. In the files, the description of a sequence is defined as `<family>|<virus name>|<protein function>`.

## Real metatranscriptomes

`real.21NBJANKW.fna`, `real.21NBJUNKW.fna`, and `real.21NBJUNYW.fna` contain the segmented RNA viral contigs identified from the Real metatranscriptomes. In the files, the description of a sequence is defined as `<family>|<length>|<method>|<cds region>|<evalue>|<virus name>|<protein function>`.