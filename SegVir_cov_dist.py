from Bio import SeqIO
import numpy as np
import os
import subprocess
import util


def run_coverm(read_path, ref_contig, threads, out_cov_path):
    """Calcualte the coverage using CoverM.
    """
    if len(read_path)==2:
        read_1 = read_path[0]
        read_2 = read_path[1]
        subprocess.run(f"coverm contig --coupled {read_1} {read_2} --reference {ref_contig} -t {threads} -o {out_cov_path}", 
                    shell=True)

    else:
        subprocess.run(f"coverm contig --single {read_path[0]} --reference {ref_contig} -t {threads} -o {out_cov_path}", 
                    shell=True)


# def cal_coverage():
#     """Calculate the coverage of sequences giving sequencing samples and reference contigs.
#     """
#     pass


def cal_pcc(vec_s1, vec_s2):
    """Calcualte the Pearson correlation coefficient (PCC) of the sets of sequences.
    """
    nor_s1 = np.array(vec_s1)/ max(vec_s1)
    nor_s2 = np.array(vec_s2)/ max(vec_s2)
    
    # pear_coff = pearsonr(nor_s1, nor_s2)[0]
    # if pear_coff < 0.0:
    #     pear_coff = 0.0
    # pear_coff = np.linalg.norm((nor_s1 - nor_s2), ord=1)/ nor_s1.shape[0]

    # pear_coff = np.linalg.norm((nor_s1 - nor_s2))/ np.sqrt(nor_s1.shape[0])
    # pear_coff = np.sum((nor_s1 - nor_s2)**2)/ nor_s1.shape[0]
    pear_coff = np.linalg.norm((nor_s1 - nor_s2))
    pear_coff = 1 - pear_coff
    if pear_coff < 0.0:
        pear_coff = 0.0

    return pear_coff


def output_coff(query_path, rdrp_path, reads_path_file, temp_dir, 
                pcc_out_path, num_thread=8):
    """Calculate the coverage.
    """
    # load all the sequences
    query_set, rdrp_set = set([]), set([])
    all_seqs = []
    for s in SeqIO.parse(query_path, 'fasta'):
        query_set.add(s.id)
        all_seqs.append(s)
    for s in SeqIO.parse(rdrp_path, 'fasta'):
        rdrp_set.add(s.id)
        all_seqs.append(s)
    SeqIO.write(all_seqs, f"{temp_dir}/all_contig.fna", 'fasta')

    # calculate coverage distribution for all the contigs
    num_samples = 0
    with open(reads_path_file) as reads_p:
        for read_index, l in enumerate(reads_p):
            reads_path = l.strip().split(',')
            run_coverm(read_path=reads_path, 
                       ref_contig=f"{temp_dir}/all_contig.fna", 
                       threads=num_thread, 
                       out_cov_path=f"{temp_dir}/{read_index}.cov.csv")
            num_samples += 1

    # load coverage distribution for all the contigs
    query_cov_dict = {s: [] for s in query_set}
    rdrp_cov_dict = {s: [] for s in rdrp_set}
    for read_index in range(num_samples):
        with open(f"{temp_dir}/{read_index}.cov.csv") as cov_rst:
            for l in cov_rst:
                l = l.strip().split('\t')
                s_id = l[0]
                if s_id != 'Contig':
                    cov = float(l[1])
                    if s_id in query_cov_dict:
                        query_cov_dict[s_id].append(cov)
                    elif s_id in rdrp_cov_dict:
                        rdrp_cov_dict[s_id].append(cov)

    # calculate the Pearson correlation coefficient (PCC) between pairs of contigs
    with open(pcc_out_path, 'w') as pcc_out:
        pcc_out.write(f"non_rdrp\trdrp\tpcc\n")
        for query_id, query_csv in query_cov_dict.items():
            temp_dict = {}
            for rdrp_id, rdrp_csv in rdrp_cov_dict.items():
                pcc = cal_pcc(query_csv, rdrp_csv)
                temp_dict[rdrp_id] = pcc
            temp_dict = {k: v for k, v in sorted(temp_dict.items(), key=lambda item: item[1], reverse=True)}
        for rdrp_id, pcc in temp_dict.items():
            pcc_out.write(f"{query_id}\t{rdrp_id}\t{pcc}\n")


if __name__ == "__main__":

    segvir_args = util.cmd_segvir_cov()
    
    rst_dir = segvir_args.rst_dir
    reads_path_file = segvir_args.reads_path
    temp_dir = segvir_args.temp_dir
    pcc_out_path = segvir_args.out_rst

    family_list = util.seg_family()
    for family in family_list:
        if os.path.exist(f"{rst_dir}/{family}.rdrp.fna"):
            if os.path.exist(f"{rst_dir}/{family}.nonrdrp.fna"):
                output_coff(query_path=f"{rst_dir}/{family}.nonrdrp.fna", 
                            rdrp_path=f"{rst_dir}/{family}.rdrp.fna", 
                            reads_path_file=reads_path_file, temp_dir=temp_dir, 
                            pcc_out_path=pcc_out_path, num_thread=8)
