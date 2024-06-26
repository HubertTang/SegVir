import argparse
import warnings
warnings.simplefilter(action='ignore', category=DeprecationWarning)


SRGVIR_VERSION = "1.0"


def seg_family():
    family_list = [ 'Alphatetraviridae', 'Alternaviridae', 'Amnoonviridae', 'Arenaviridae', 'Aspiviridae', 
                    'Benyviridae', 'Birnaviridae', 'Bromoviridae', 'Chrysoviridae', 'Chuviridae', 
                    'Closteroviridae', 'Cruliviridae', 'Cystoviridae', 'Discoviridae', 'Fimoviridae', 
                    'Hadakaviridae', 'Hantaviridae', 'Kitaviridae', 'Leishbuviridae', 'Mayoviridae', 
                    'Megabirnaviridae', 'Nairoviridae', 'Nodaviridae', 'Orthomyxoviridae', 'Partitiviridae', 
                    'Peribunyaviridae', 'Phasmaviridae', 'Phenuiviridae', 'Picobirnaviridae', 'Polymycoviridae', 
                    'Potyviridae', 'Qinviridae', 'Quadriviridae', 'Rhabdoviridae', 'Secoviridae', 
                    'Sedoreoviridae', 'Spinareoviridae', 'Tospoviridae', 'Virgaviridae', 'Yueviridae']
    return family_list


def cmd_segvir():
    parser = argparse.ArgumentParser(description="SegVir: identify \033[4mSeg\033[0mmented rna \033[4mVir\033[0muses and reconstruct their complete genomes from metatranscriptomes.")

    # argument for dataset
    parser.add_argument(
        '--input',
        type=str,
        required=True, 
        help="Path of the query contigs (in 'fasta' format)."
    )

    parser.add_argument(
        '--outdir',
        type=str,
        required=True, 
        help="Directory to store results. The directory will be created if it does not exist."
    )

    parser.add_argument(
        "--database",
        type=str,
        default="segvir_db",
        help="The database directory. (Use the absolute path to specify the location of the database. Default: SegVir/segvir_db)"
        )

    parser.add_argument(
        "--min_len",
        type=int,
        default=300,
        help="The minimal length of the contigs (default: 300nt)."
        )

    parser.add_argument(
        "--host",
        default=None,
        type=str,
        help="Path of the host genome. (The path can be a fasta file or a directory containing the genomes.)")
    
    parser.add_argument(
        "--blastp",
        default=1e-5,
        type=float,
        help="The minimun e-value of BLASTP (default: 1e-5).")
    
    parser.add_argument(
        "--hmmer",
        default=1e-5,
        type=float,
        help="The minimun e-value of HMMER  (default: 1e-5).")
    
    parser.add_argument(
        "--rdrp_evalue",
        default=1e-5,
        type=float,
        help="The minimun e-value of identified RdRp (default: 1e-10).")
    
    parser.add_argument(
        "--rdrp_len",
        default=900,
        type=float,
        help="The minimun length of identified RdRp  (default: 900bp).")

    parser.add_argument(
        "--vote_thres",
        default=0.85,
        type=float,
        help="Multiply the best hit's bitscore by this parameter as the threshold, and determine the taxonomy by majority voting on the results above the threshold.  (default: 0.85).")

    parser.add_argument(
        '-t', "--thread",
        default=8,
        type=int,
        help="The number of threads  (default: 8).")

    parser.add_argument(
        "--tempdir",
        type=str,
        default='temp',
        help="The temporary directory (default: <outdir>/temp)."
        )

    parser.add_argument(
        "--outfmt",
        type=int,
        default=1,
        help="The output format of identified viral genomes (default: 1).\
            1: save all the sequences in a file.\
            2: save the sequences by families.\
            3. save the sequences by families and RdRp."
        )

    # version
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'SegVir_v{SRGVIR_VERSION}'
    )

    args = parser.parse_args()

    return args


def cmd_segvir_polish():
    parser = argparse.ArgumentParser(description="SegVir: identify \033[4mSeg\033[0mmented rna \033[4mVir\033[0muses and reconstruct their complete genomes from metatranscriptomes.")

    parser.add_argument(
        "--tempdir",
        type=str,
        required=True, 
        help="The temporary directory generated by running `SegVir.py`."
        )

    parser.add_argument(
        '--outdir',
        type=str,
        required=True, 
        help="Directory to store results. If exists, the existing results in directory will be overwrite if set `--overwrite` as True.\
             If not exist, the directory will be created."
    )

    parser.add_argument(
        "--database",
        type=str,
        default="segvir_db",
        help="The database directory. (Use the absolute path to specify the location of the database. Default: SegVir/segvir_db)"
        )

    parser.add_argument(
        "--min_len",
        type=int,
        default=300,
        help="The minimal length of the contigs (default: 300nt)."
        )

    parser.add_argument(
        "--blastp",
        default=1e-5,
        type=float,
        help="The minimun e-value of BLASTP (default: 1e-5).")
    
    parser.add_argument(
        "--hmmer",
        default=1e-5,
        type=float,
        help="The minimun e-value of HMMER  (default: 1e-5).")
    
    parser.add_argument(
        "--rdrp_evalue",
        default=1e-5,
        type=float,
        help="The minimun e-value of identified RdRp (default: 1e-10).")
    
    parser.add_argument(
        "--rdrp_len",
        default=900,
        type=float,
        help="The minimun length of identified RdRp  (default: 900bp).")

    parser.add_argument(
        "--vote_thres",
        default=0.85,
        type=float,
        help="Multiply the best hit's bitscore by this parameter as the threshold, and determine the taxonomy by majority voting on the results above the threshold.  (default: 0.85).")

    parser.add_argument(
        "--overwrite",
        type=bool,
        default=False,
        help="Overwrite the existing results in the output directory (default: False)."
        )

    parser.add_argument(
        "--outfmt",
        type=int,
        default=1,
        help="The output format of identified viral genomes (default: 1).\
            1: save all the sequences in a file.\
            2: save the sequences by families.\
            3. save the sequences by families and RdRp."
        )
    
    # version
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'SegVir_v{SRGVIR_VERSION}'
    )

    args = parser.parse_args()

    assert args.outfmt in [1, 2, 3], 'Unknown output format (seq_outfmt: 1 (default), 2, 3).'

    return args


def cmd_segvir_cov():
    parser = argparse.ArgumentParser(description="SegVir: identify \033[4mSeg\033[0mmented rna \033[4mVir\033[0muses and reconstruct their complete genomes from metatranscriptomes.")

    parser.add_argument(
        '--rst_dir',
        type=str,
        required=True, 
        help="Path of the output RdRp and non-RdRp contigs (in 'fasta' format)."
    )

    parser.add_argument(
        '--reads_path',
        type=str,
        required=True, 
        help="File storing the path of metagenomic samples."
    )

    parser.add_argument(
        '--out_rst',
        type=str,
        required=True, 
        help="The path of output results."
    )

    parser.add_argument(
        '--temp_dir',
        type=str,
        required=True, 
        help="The temporary directory."
    )

    # version
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'SegVir_v{SRGVIR_VERSION}'
    )

    args = parser.parse_args()

    return args


