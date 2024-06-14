# SegVir
SegVir is designed to identify and reconstruct complete genomes of segmented RNA viruses from complex metatranscriptomic data. The tool utilizes both close and remote homology searches to detect conserved and divergent viral segments, and introduces a new method to evaluate genome completeness based on protein clusters with similar functions.

## Required Dependencies

* Python 3.x
* diamond
* blast
* hmmer
* biopython
* pandas
* ORFfinder

## Quick install (Linux only)

1. Download SegVir by `git clone`

   ```bash
   git clone https://github.com/HubertTang/SegVir.git
   cd SegVir
   ```

2. We recommend using `conda` to install all the dependencies.

   ```bash
   # install the plasme
   conda env create -f segvir.yaml
   # install ORFfinder
   conda activate segvir
   bash install_orffinder.sh
   conda deactivate
   # activate the environment
   conda activate segvir
   ```

3. Download the reference database

   Download the reference database (226MB) from [OneDrive](https://portland-my.sharepoint.com/:u:/g/personal/xubotang2-c_my_cityu_edu_hk/EYRIkHnE58xIrWH1tBzPl_MBK0DNx4YfIf8IVhpmwUzk4g?e=iTSiDY) to the same directory with `SegVir.py` and uncompress it. If you uncompress the database somewhere else, please specify that path in `--database`.

## Usage

SegVir requires input assembled contigs in fasta format, and outputs include the sequences and detailed information of the identified segmented RNA virus genomes.

```bash
python SegVir.py --input [INPUT_CONTIG] --outdir [OUTPUT_DIRECTORY] [OPTIONS]
```

Required arguments: 

   `--input`: Path of the query contigs (in 'fasta' format).

   `--outdir`: Directory to store results. The directory will be created if it does not exist.

 More optional arguments:

   `--database`: The database directory. (Use the absolute path to specify the location of the database. Default: SegVir/segvir_db)

   `--min_len`: The minimal length of the contigs (default: 300bp).

   `--host`: Path of the host genome. (The path can be a fasta file or a directory containing the genomes.)

   `--blastp`: The minimun e-value of BLASTP (default: 1e-5).

   `--hmmer`: The minimun e-value of HMMER (default: 1e-5).

   `--rdrp_evalue`: The minimun e-value of identified RdRp (default: 1e-10).

   `--rdrp_len`: The minimun length of identified RdRp (default: 900bp).

   `--tempdir`: The temporary directory (default: \<outdir>/temp).

   `--vote_thres`: Multiply the best hit's bitscore by this parameter as the threshold, and determine the taxonomy by majority voting on the results above the threshold.  (default: 0.85).

   `--outfmt`: The output format of identified viral genomes (default: 1). 1: save the identified genomes in a file; 2: save the identified genomes by families; 3. save the identified genomes by families and RdRp."

   `-t`, `--thread`: The number of threads (default: 8).

## Outputs

### Output files

| Files            | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| segvir.fna       | Fasta file of all identified segmented RNA virus genomes     |
| segvir.csv       | Report file of the details of the identified segmented RNA virus genomes |
| segvir.score.csv | Report file of the completeness and conservation score of the identified viral genomes |

### Output report format

| `segvir.csv` | Description                                                  |
| ------------ | ------------------------------------------------------------ |
| family       | The taxon of identified contig in the family level           |
| contig       | The ID of the contig                                         |
| length       | The length of the contig                                     |
| gene_coor    | The coordinate of the identified gene                        |
| blastp_e     | The e-value outputted by Diamond BLASTp                      |
| hmmer_e      | The e-value outputted by HMMER                               |
| ref_name     | The organism name of the aligned reference sequence using BLASTp. |
| function     | The function of the aligned reference sequence/ function cluster |
| fc           | The ID of the aligned function cluster                       |

| `segvir.score.csv` | Description                                                  |
| ------------------ | ------------------------------------------------------------ |
| family             | Family of the identified segmented RNA genomes               |
| completeness       | The estimated completeness of the identified viral genome    |
| cs                 | The conservation score of the identified viral genome        |
| ref_cs             | The reference conservation score of the closet reference genome |

## Example

```bash
python SegVir.py --input example/test.fna --outdir example/out --host example/ref_hosts.fna
```

Identify the segmented RNA viruses from  `test.fna` in `example` directory. The host genomes is saved in `example/ref_hosts.fna` and the results will be saved in `example/out`. 

## (Optional) Polish the results

You may not be satisfied with the results you get the first time, for example, the set thresholds are too loose and the results contain many potential false positives. You can run `SegVir_polish.py` to reset the thresholds and output the corresponding results without re-running `SegVir.py` (alignment is very time consuming).

```bash
python SegVir_polish.py --tempdir [TEMP_DIR] --outdir [OUTPUT_DIRECTORY] [OPTIONS]
```

Required arguments: 

   `--tempdir`: The temporary directory generated by running `SegVir.py`.

   `--outdir`: Directory to store results. If exists, the existing results in directory will be overwrite if set `--overwrite` as True. If not exist, the directory will be created.

 More optional arguments:

   `--database`: The database directory. (Use the absolute path to specify the location of the database. Default: SegVir/segvir_db)

   `--min_len`: The minimal length of the contigs (default: 300bp).

   `--blastp`: The minimun e-value of BLASTP (default: 1e-5).

   `--hmmer`: The minimun e-value of HMMER (default: 1e-5).

   `--rdrp_evalue`: The minimun e-value of identified RdRp (default: 1e-10).

   `--rdrp_len`: The minimun length of identified RdRp (default: 900bp).

   `--vote_thres`: Multiply the best hit's bitscore by this parameter as the threshold, and determine the taxonomy by majority voting on the results above the threshold.  (default: 0.85).

   `--outfmt`: The output format of identified viral genomes (default: 1). 1: save the identified genomes in a file; 2: save the identified genomes by families; 3. save the identified genomes by families and RdRp."

   `--overwrite`: Overwrite the existing results in the output directory (default: False).

**Example:**

```bash
python SegVir_polish.py --tempdir example/out/temp --outdir example/out_polish --rdrp_evalue 1e-100 --outfmt 3 --overwrite True
```

Polish the results of the first use of using `SegVir.py` to identify segmented RNA viral genomes in  `example/test.fna`. The polished results will be saved in `example/out_polish`.

## (Optional) Calculate the Pearson correlation coefficient of coverage distributions between non-RdRp and RdRp contigs

The script `SegVir_cov_dist.py` can calculate the correlation (Pearson correlation coefficient) between non-RdRp contigs and RdRp contigs by using the coverage distribution of contigs across multiple samples. To do this, you need to install `coverm` using the following command:

```bash
conda install coverm
```

The usage of the script is:

```bash
python SegVir_cov_dist.py --rst_dir [RST_DIR] --reads_path [READS_PATH] --out_rst [OUT_RST] --temp_dir [TEMP_DIR]
```

Required arguments: 

`--rst_dir`: the path of the directory for the results generated using `--outfmt 3`.

`--reads_path`: file storing the path of sequencing samples.

If it is single-end reads, the file format is as follows:

```
read/path/1
read/path/2
... ...
```

If it is pair-end reads, the file format is as follows:

```
read/path/1_1,read/path/1_2
read/path/2_1,read/path/2_2
... ...
```

`-out_rst`: the path of output results.

`--temp_dir`: the temporary directory.

## Supporting data

The `supp_data` directory saves the identified segmented RNA viral contigs from `Missing segment datasets` and `Real metatranscriptomes`.
