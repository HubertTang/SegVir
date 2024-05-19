from Bio import SeqIO
import math
from multiprocessing import Pool
import os
import pandas as pd
import random
import shutil
import util
import subprocess


FAMILY_LIST = util.seg_family()


def build_host_db(host_path, temp_dir):
    """Check if the database is complete and build the host database.
    """
    # pre-build the host database
    if os.path.isfile(host_path):
        subprocess.run(f"makeblastdb -in {host_path} -dbtype nucl -out {temp_dir}/host",
                    shell=True)
    else:
        with os.scandir(host_path) as it:
            if any(it):
                print('The directory of host genome is empty.')
                exit(0)

        with open(f"{temp_dir}/host.fna",'wb') as wfd:
            for f in os.listdir(host_path):
                with open(f"{host_path}/{f}",'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
        subprocess.run(f"makeblastdb -in {temp_dir}/host.fna -dbtype nucl -out {temp_dir}/host",
                    shell=True)
        os.remove(f"{temp_dir}/host.fna")


def split_fasta(fasta_file, num_split=10):
    """Split original fasta file into several fasta files.
    """
    all_seq_id = [s.id for s in SeqIO.parse(fasta_file, 'fasta')]
    random.shuffle(all_seq_id)

    num_per_split = math.ceil(len(all_seq_id)/ num_split)

    num_files = 0
    seq_index = SeqIO.index(fasta_file, 'fasta')
    for i in range(0, len(all_seq_id), num_per_split):
        if len(all_seq_id[i: i + num_per_split]) > 0:
            out_seq = [seq_index[s] for s in all_seq_id[i: i + num_per_split]]
            SeqIO.write(out_seq, f"{fasta_file}.{num_files}", 'fasta')
            num_files += 1

    return num_files


def orffinder(dna_path):
    """Run ORFfinder.
    """
    subprocess.run(f"./ORFfinder -n True -in {dna_path} -out {dna_path}.aa", shell=True)


def run_multi_orffinder(contig_path, threads=32):
    """Run ORFfinder in multiple threads.
    """
    num_files = split_fasta(fasta_file=contig_path, num_split=threads)
    threads = num_files
    pool = Pool(processes=threads)
    for temp_id in range(threads):
        pool.apply_async(orffinder, [f"{contig_path}.{temp_id}"])
    pool.close()
    pool.join()
    
    # need to merge and delete the temporary files
    # work_dir = os.path.dirname(contig_path)
    for i in range(threads):
        os.remove(f"{contig_path}.{i}")
    
    # merge the files
    with open(f'{contig_path}.aa', 'w') as outfile:
        for i in range(threads):
            fname = f"{contig_path}.{i}.aa"
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            os.remove(fname)


def run_diamond(contig_path, db_path, diamond_out, temp_dir, num_thread=56):
    """Run DIAMOND BLASTx.
    """
    subprocess.run(f"diamond blastp --threads {num_thread} -d {db_path} -q {contig_path} -o {diamond_out} --sensitive -c1 --tmpdir {temp_dir}", shell=True)


def parse_diamond(diamond_out_path, e_thres=0.001):
    """Parse DIAMOND BLASTx output.
    """
    # parse the diamond results, use the best hit
    # output result dictionary {nt: {prot: [subject, evalue, prot_start, prot_end]}}
    query_dict = {}
    with open(diamond_out_path) as bo:
        for l in bo:
            records = l.split()
            prot_id = records[0]
            nt_id = records[0].split('_', 1)[1].split(':')[0]
            ref_id = records[1]
            e_value = float(records[10])
            q_start = int(records[6])
            q_end = int(records[7])
            
            if e_value <= e_thres:
                if nt_id not in query_dict:
                    query_dict[nt_id] = {}
                    query_dict[nt_id][prot_id] = [ref_id, e_value, q_start, q_end]
                else:
                    if prot_id in query_dict[nt_id]:
                        if e_value < query_dict[nt_id][prot_id][1]:
                            query_dict[nt_id][prot_id] = [ref_id, e_value, q_start, q_end]
                    else:
                        query_dict[nt_id][prot_id] = [ref_id, e_value, q_start, q_end]

    return query_dict


def run_hmm(hmm_path, prot_path, out_path, num_thread=56):
    """Run HMM on the target proteins.
    """
    # subprocess.run(f"hmmsearch --noali --cut_nc -o {out_path}.hmm --domtblout {out_path}.tblout --cpu {num_thread} {hmm_path} {prot_path}", shell=True)
    subprocess.run(f"hmmsearch --noali -o {out_path}.hmm --domtblout {out_path}.tblout --cpu {num_thread} {hmm_path} {prot_path}", shell=True)


def parse_hmmer(hmm_out_path):
    """Parse HMMER output.
    """
    # parse the HMMER result, use the best hit, 
    # output result dictionary {nt: {prot: [pc, evalue, prot_start, prot_end]}}
    query_dict = {}
    with open(hmm_out_path) as hop:
        for l in hop:
            if l[0] != '#' and l[0] != '':
                l = l.strip().split()
                prot_id = l[0]
                nt_id = prot_id.split('_', 1)[1].split(':')[0]
                ref_id = l[3]
                e_value = float(l[6])
                q_start = int(l[17])
                q_end = int(l[18])
                
                if nt_id not in query_dict:
                    query_dict[nt_id] = {}
                    query_dict[nt_id][prot_id] = [ref_id, e_value, q_start, q_end]
                else:
                    if prot_id in query_dict[nt_id]:
                        if e_value <= query_dict[nt_id][prot_id][1]:
                            query_dict[nt_id][prot_id] = [ref_id, e_value, q_start, q_end]
                    else:
                        query_dict[nt_id][prot_id] = [ref_id, e_value, q_start, q_end]
        
    return query_dict


def merge_range(range_list):
    saved = list(range_list[0])
    for st, en in sorted([sorted(t) for t in range_list]):
        if st <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)


def parse_remove_host(contig_path, contig_aa_path, diamond_dict, hmm_nucl_dict, temp_dir):
    """Remove the host contigs, and identify contigs containing genes from hosts and virus.
    Due to the presence of nrEVEs and chimeric assembly contigs, it is necessary to identify and detect whether the detected virus may originate from nrEVE.
    diamond dictionary {query: [fam, segment, function, evalue, query start, query end, reference]}
    hmmer dictionary (query: [fam, segment, function, evalue, query start, query end, query_trans_id])
    """
    # build the index of contigs
    seq_index = SeqIO.index(contig_path, 'fasta')
    
    # save the identified contigs (DIAMOND & HMMER) in the temporary directory
    indent_seq_set = set()
    for s in diamond_dict:
        indent_seq_set.add(s)
    for s in hmm_nucl_dict:
        indent_seq_set.add(s)
    out_seq_list = []
    for s in indent_seq_set:
        out_seq_list.append(seq_index[s])
    SeqIO.write(out_seq_list, f"{temp_dir}/temp_seq.fna", 'fasta')

    # run BLASTN to remove the host contigs
    with open(f"{temp_dir}/contig.blastn", 'w') as blastn_out:
        with open(f"{temp_dir}/contig.host.blastn") as blastn_in:
            for l in blastn_in:
                query = l.split('\t')[0]
                if query in indent_seq_set:
                    blastn_out.write(l)

    # build contig to length dictionary
    contig_len_dict = {}
    for s in SeqIO.parse(contig_path, 'fasta'):
        contig_len_dict[s.id] = len(s.seq)

    # build contig to translated dictionary
    contig_prot_coor_dict = {}
    for s in SeqIO.parse(contig_aa_path, 'fasta'):
        trans_start = int(s.id.split(':')[-2])
        trans_end = int(s.id.split(':')[-1])
        contig_prot_coor_dict[s.id] = (trans_start - 1, trans_end)

    # analyze the blastn result (the best hit)
    # need to count all the alignment region for the host in order to distinguish the chimeric assembly contigs and nrEVE
    contig_aln_dict = {}
    with open(f"{temp_dir}/contig.blastn") as blas_p:
        for l in blas_p:
            records = l.split('\t')
            query = records[0]
            reference = records[1]
            ident = float(records[2])
            qstart = int(records[6])
            qend = int(records[7])
            evalue = float(records[10])
            overall_cov = int(records[12])
            if query not in contig_aln_dict:
                contig_aln_dict[query] = [ident, qstart, qend, evalue, overall_cov, reference, [(qstart, qend)]]
            else:
                if reference == contig_aln_dict[query][5]:
                    # if ident >= 80 and evalue <= 1e-5:
                    contig_aln_dict[query][6] = list(merge_range(contig_aln_dict[query][6] + [(qstart, qend)]))

    # analyze if the contigs are from host or nrEVEs/ chimeric assembly contigs
    # Strategy for determining the origins of contig s:
    # 1) host: a) s is aligned to hosts, not aligned to viruses 
    #          b) s is aligned to hosts (e-value <= 1e-5) and aligned to viruses (e-value >= 1e-5)
    # 2) nrEVEs/ chimeric: s is aligned to both hosts and viruses (e-value <= 1e-5)

    # output the alignment results 
    with open(f"{temp_dir}/candidate_blast.host.csv", 'w') as host_out:
        host_out.write(f"query\treference\taln_start\taln_end\n")
        for query, (ident, qstart, qend, evalue, overall_cov, reference, aln_coor_list) in contig_aln_dict.items():
            for (q_start, q_end) in aln_coor_list:
                host_out.write(f"{query}\t{reference}\t{q_start}\t{q_end}\n")


def tuple_list2set(tuple_list):
    """Convert the list of tuples to a set of numbers.
    """
    out_set = set([])
    for (t1, t2) in tuple_list:
        for i in range(t1, t2+1):
            out_set.add(i)
    return out_set


def run_segvir_aln(contig_path, host_db_path, ref_dir, temp_dir, threads=8, len_thres=300):
    """Run SegVir on the simulated sequencing datasets.
    """
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    contig_aa_path = f"{temp_dir}/contig.clean.fna.aa"
    contig_seq_index = {}
    for s in SeqIO.parse(contig_path, 'fasta'):
        contig_seq_index[s.id] = s
    # contig_aa_index = SeqIO.index(contig_aa_path, 'fasta')
    contig_len_dict = {s: len(contig_seq_index[s]) for s in contig_seq_index}

    ##### Step 0: run BLASTN to remove host contigs (95% identity and 95% coverage) and predict the proteins from contigs
    with open(f"{temp_dir}/contig.fna", 'w') as filter_out:
        for s, s_len in contig_len_dict.items():
            if s_len >= len_thres:
                filter_out.write(f">{s}\n{contig_seq_index[s].seq}\n")

    if os.stat(f"{temp_dir}/contig.fna").st_size == 0:
        print(f"No sequence is longer than {len_thres}.")
        exit(0)

    if host_db_path != None:
        subprocess.run(f"blastn -db {host_db_path} -query {temp_dir}/contig.fna -num_threads {threads} -out {temp_dir}/contig.host.blastn \
                        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp'", shell=True)

        host_contig_set = set([])
        with open(f"{temp_dir}/contig.host.csv", 'w') as host_out:
            host_out.write(f"query\treference\tidentity\tcoverage\n")
            with open(f"{temp_dir}/contig.host.blastn") as blastn_out:
                for l in blastn_out:
                    l = l.strip().split()
                    query = l[0]
                    if query not in host_contig_set:
                        sub = l[1]
                        ident = float(l[2])
                        cov = float(l[12])
                        if ident >= 95 and cov >= 95:
                            host_out.write(f"{query}\t{sub}\t{ident}\t{cov}\n")
                            host_contig_set.add(query)

        with open(f"{temp_dir}/contig.clean.fna", 'w') as clean_out:
            for s in SeqIO.parse(f"{temp_dir}/contig.fna", 'fasta'):
                if s.id not in host_contig_set:
                    clean_out.write(f">{s.description}\n{s.seq}\n")
    else:
        shutil.copy(f"{temp_dir}/contig.fna", f"{temp_dir}/contig.clean.fna")
        open(f"{temp_dir}/contig.host.blastn", 'a').close()

    print(f"Predict the proteins using ORFfinder from the input contigs ... ...")
    # subprocess.run(f"./ORFfinder -in {temp_dir}/contig.clean.fna -out {contig_aa_path}", shell=True)
    run_multi_orffinder(contig_path=f"{temp_dir}/contig.clean.fna", threads=threads)

    # # Step 1: load the taxonomy information
    # # Load the taxon information
    # print(f"Load the taxonomy ... ...")
    # taxon_dir = f"{ref_dir}/taxon"
    # acc_name_dict = {}
    # for f in FAMILY_LIST:
    #     input_csv = f"{taxon_dir}/{f}.csv"
    #     with open(input_csv) as csv_in:
    #         for l in csv_in:
    #             l_info = l.strip().split('\t')
    #             s_id, taxon_name = l_info[0], l_info[1]
    #             acc_name_dict[s_id] = [f, taxon_name]
        
    ##### Step 1: run DIAMOND BLASTp
    print(f"Run diamond ... ...")
    run_diamond(contig_path=contig_aa_path, db_path=f"{ref_dir}/ref", 
                diamond_out=f"{temp_dir}/contig.diamond", 
                temp_dir=temp_dir, num_thread=threads)

    # Step 2: run HMMER
    print(f"Run HMMER ... ...")
    run_hmm(hmm_path=f"{ref_dir}/seg_non_seg.temp.hmm",
            prot_path=contig_aa_path,
            out_path=f"{temp_dir}/hmm",
            num_thread=threads)
    
    # # if `consise_rst=True`, keep the contigs with the smallest e-value
    # if consise_rst:
    #     pass
    # else:
    #     os.rename(f'{temp_dir}/contig.temp.diamond', f'{temp_dir}/contig.diamond')
    #     os.rename(f'{temp_dir}/hmm.temp.tblout', f'{temp_dir}/hmm.temp.tblout')


def segvir_parse_rst(contig_path, ref_dir, temp_dir, blastp_evalue=1e-5, hmmer_evalue=1e-5, host_ratio_thres=0.84,
                     rdrp_evalue=1e-10, rdrp_len=900):
    """Parse the alignment and output the final result.
    """
    # load the translated proteins
    contig_aa_path = f"{temp_dir}/contig.clean.fna.aa"
    
    # load the index of all contigs
    contig_seq_index = {}
    for s in SeqIO.parse(contig_path, 'fasta'):
        contig_seq_index[s.id] = s
    contig_len_dict = {s: len(contig_seq_index[s]) for s in contig_seq_index}

    # load the Diamond and HMMER results
    diamond_dict = parse_diamond(diamond_out_path=f"{temp_dir}/contig.diamond")
    hmm_dict = parse_hmmer(hmm_out_path=f"{temp_dir}/hmm.tblout")
    
    # load the taxonomy information
    taxon_dir = f"{ref_dir}/taxon"
    acc_name_dict = {}
    for f in FAMILY_LIST:
        input_csv = f"{taxon_dir}/{f}.csv"
        with open(input_csv) as csv_in:
            for l in csv_in:
                l_info = l.strip().split('\t')
                s_id, taxon_name = l_info[0], l_info[1]
                acc_name_dict[s_id] = [f, taxon_name]

    # analyze if the contigs are from host or nrEVEs/ chimeric assembly contigs
    # Strategy for determining the origins of contig s:
    # 1) host: a) s is aligned to hosts, not aligned to viruses 
    #          b) s is aligned to hosts (e-value <= 1e-5) and aligned to viruses (e-value >= 1e-5)
    # 2) nrEVEs/ chimeric: s is aligned to both hosts and viruses (e-value <= 1e-5)
    # load the aligned region of contigs and calculate the aligned ratio with hosts
    diamond_coor_dict, hmmer_coor_dict, host_coor_dict = {}, {}, {}
    # {nt: {prot: [ref, evalue, prot_start, prot_end]}}
    with open(f"{temp_dir}/candidate_diamond.virus.csv", 'w') as diamond_out:
        diamond_out.write(f"query\tprotein\treference\tevalue\torf_start\torf_end\n")
        for query, prot_dict in diamond_dict.items():
            for prot_id, (ref, evalue, prot_start, prot_end) in prot_dict.items():
                orf_start = min(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
                orf_end = max(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
                diamond_out.write(f"{query}\t{prot_id}\t{ref}\t{evalue}\t{orf_start}\t{orf_end}\n")
    
    diamond_df = pd.read_csv(f"{temp_dir}/candidate_diamond.virus.csv", sep='\t')
    for query, ref, orf_s, orf_e in zip(diamond_df['query'], diamond_df['reference'], diamond_df['orf_start'], diamond_df['orf_end']):
        if query not in diamond_coor_dict:
            diamond_coor_dict[query] = [(orf_s, orf_e)]
        else:
            diamond_coor_dict[query] = list(merge_range(diamond_coor_dict[query] + [(orf_s, orf_e)]))

    # {nt: {prot: [pc, evalue, prot_start, prot_end]}}
    with open(f"{temp_dir}/candidate_hmm.virus.csv", 'w') as hmmer_out:
        hmmer_out.write(f"query\tprotein\treference\tevalue\torf_start\torf_end\n")
        for query, prot_dict in hmm_dict.items():
            for prot_id, (ref, evalue, prot_start, prot_end) in prot_dict.items():
                orf_start = min(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
                orf_end = max(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
                hmmer_out.write(f"{query}\t{prot_id}\t{ref}\t{evalue}\t{orf_start}\t{orf_end}\n")

    hmmer_df = pd.read_csv(f"{temp_dir}/candidate_hmm.virus.csv", sep='\t')
    for query, ref, orf_s, orf_e in zip(hmmer_df['query'], hmmer_df['reference'], hmmer_df['orf_start'], hmmer_df['orf_end']):
        if ref.split('_')[0] != 'FP':
            if query not in hmmer_coor_dict:
                hmmer_coor_dict[query] = [(orf_s, orf_e)]
            else:
                hmmer_coor_dict[query] = list(merge_range(hmmer_coor_dict[query] + [(orf_s, orf_e)]))

    ##### Step 3: remove/ detect the host contamination
    # merge the results of DIAMOND, HMM, and CAT to obtain 'temp_seq.virus.potential.fna'
    print(f"Remove the host contamination ... ...")  
    parse_remove_host(contig_path=f"{temp_dir}/contig.clean.fna", contig_aa_path=contig_aa_path,
                        diamond_dict=diamond_dict, hmm_nucl_dict=hmm_dict,
                        temp_dir=temp_dir)
    host_df = pd.read_csv(f"{temp_dir}/candidate_blast.host.csv", sep='\t')
    for query, ref, aln_s, aln_e in zip(host_df['query'], host_df['reference'], host_df['aln_start'], host_df['aln_end']):
        if query not in host_coor_dict:
            host_coor_dict[query] = [(aln_s, aln_e)]
        else:
            host_coor_dict[query] = list(merge_range(host_coor_dict[query] + [(aln_s, aln_e)]))

    # Step 4: integrate the results of all the candidate viral contigs (output diamond first, then hmmer)
    with open(f"{temp_dir}/candidate_result.csv", 'w') as cand_out:
        cand_out.write(f"query\tprot\tmethod\treference\tevalue\thost_ratio\n")
        for prot in SeqIO.index(contig_aa_path, 'fasta'):
            query = prot.split('_', 1)[1].split(':')[0]
            
            # count the aligned regions based on hosts
            host_coor_set = set([])
            if query in host_coor_dict:
                host_coor_set = tuple_list2set(host_coor_dict[query])

            if query in diamond_dict:
                diamond_coor_set = tuple_list2set(diamond_coor_dict[query])
                diamond_overlap_ratio = len(host_coor_set.intersection(diamond_coor_set))/ len(diamond_coor_set)
                if prot in diamond_dict[query]:
                    ref = diamond_dict[query][prot][0]
                    evalue = diamond_dict[query][prot][1]
                    cand_out.write(f"{query}\t{prot}\tdiamond\t{ref}\t{evalue}\t{diamond_overlap_ratio}\n")

            if query in hmmer_coor_dict:
                hmm_coor_set = tuple_list2set(hmmer_coor_dict[query])
                hmmer_overlap_ratio = len(host_coor_set.intersection(hmm_coor_set))/ len(hmm_coor_set)
                if prot in hmm_dict[query]:
                    ref = hmm_dict[query][prot][0]
                    if ref.split('_')[0] != 'FP':
                        evalue = hmm_dict[query][prot][1]
                        cand_out.write(f"{query}\t{prot}\thmmer\t{ref}\t{evalue}\t{hmmer_overlap_ratio}\n")

    # Step 5: calculate the completeness and conservation score
    # load the dictionary from PC to fc
    pc2fc_dict = {f: {} for f in FAMILY_LIST}
    fc2func_dict = {f: {} for f in FAMILY_LIST}
    fc_info_df = pd.read_csv(f"{ref_dir}/pc2fc.tsv", sep='\t')
    for family, pc, fc, anno_s, func in zip(fc_info_df['family'], fc_info_df['pc'], fc_info_df['fc'], fc_info_df['anno_state'], fc_info_df['function']):
        if pd.isna(func):
            func = 'unknown'
        pc2fc_dict[family][pc] = fc
        fc2func_dict[family][fc] = [anno_s, func]

    # load the reference fc sets and conservation score
    fc_score_dict = {f: {} for f in FAMILY_LIST}
    fc_score_df = pd.read_csv(f"{ref_dir}/ref_fc_set_cs.tsv", sep='\t')
    for fam, fc_set, score in zip(fc_score_df['family'], fc_score_df['fc_set'], fc_score_df['score']):
        fc_set = frozenset([i for i in fc_set.split(',')])
        if fc_set not in fc_score_dict[fam]:
            fc_score_dict[fam][fc_set] = score
    
    # load the fc weights
    fc2weight_dict = {f: {} for f in FAMILY_LIST}  # family: {rdrp_fc_id: {fc_id: weights}}
    weight_df = pd.read_csv(f"{ref_dir}/fc_weight.csv", sep='\t')
    for fam, rdrp_fc, fc, weight in zip(weight_df['family'], weight_df['rdrp_fc'], weight_df['fc'], weight_df['weight']):
        if rdrp_fc in fc2weight_dict[fam]:
            fc2weight_dict[fam][rdrp_fc][fc] = weight
        else:
            fc2weight_dict[fam][rdrp_fc] = {fc: weight}

    # output the results step 1: parse the results (detect if the group contains RdRp)
    # load the singleton information (protein: function)
    singleton2func_dict = {}    # singleton: {anno_s, func}
    singleton_df = pd.read_csv(f"{ref_dir}/singleton.csv", sep='\t')
    for prot_id, func in zip(singleton_df['singleton'], singleton_df['function']):
        anno_s = '-'
        if any(i in func.lower() for i in ['rdrp', 'polymerase']):
            anno_s = '*'
        singleton2func_dict[prot_id] = [anno_s, func]

    # load the single-segment viruses
    vir_single_seg_dict = {}
    single_vir_df = pd.read_csv(f"{ref_dir}/virus.single_seg.csv", sep='\t')
    for fam, vir_name in zip(single_vir_df['family'], single_vir_df['virus']):
        if fam in vir_single_seg_dict:
            vir_single_seg_dict[fam].add(vir_name)
        else:
            vir_single_seg_dict[fam] = set([vir_name])

    fam_rdrp_set = set([])
    coor2ref_dict = {}
    group_diamond_rst_dict = {f: {} for f in FAMILY_LIST}   # family: {contig: {coor: reference}}
    group_hmmer_rst_dict = {f: {} for f in FAMILY_LIST} # family: {contig: {coor: fc_index}}
    cand_rst_df = pd.read_csv(f"{temp_dir}/candidate_result.csv", sep='\t')
    for query, prot, method, ref, evalue, host_r in zip(cand_rst_df['query'], cand_rst_df['prot'], cand_rst_df['method'], cand_rst_df['reference'], cand_rst_df['evalue'], cand_rst_df['host_ratio']):
        if host_r <= host_ratio_thres:
            if method == 'diamond':
                if evalue <= blastp_evalue:
                    fam, taxa = acc_name_dict[ref][0], acc_name_dict[ref][1]
                    
                    # check if taxa is from single-segment virus
                    if fam in vir_single_seg_dict:
                        if taxa in vir_single_seg_dict[fam]:
                            if evalue <= 1e-25:
                                continue
                    
                    # check if rdrp satisfy the requirements
                    try:
                        if singleton2func_dict[ref][0] == '*':
                            if evalue > rdrp_evalue or contig_len_dict[query] < rdrp_len:
                                continue
                    except KeyError:
                        pass

                    aln_coor = prot.split(':', 1)[1]
                    if query in group_diamond_rst_dict[fam]:
                        group_diamond_rst_dict[fam][query][prot] = ref
                    else:
                        group_diamond_rst_dict[fam][query] = {prot: ref}

                    try:
                        if singleton2func_dict[ref][0] == '*':
                            if evalue <= rdrp_evalue:
                                if contig_len_dict[query] >= rdrp_len:
                                    fam_rdrp_set.add(fam)
                    except KeyError:
                        pass

                    if query in coor2ref_dict:
                        coor2ref_dict[query][prot.split(':', 1)[1]] = ref
                    else:
                        coor2ref_dict[query] = {prot.split(':', 1)[1]: ref}

            elif method == 'hmmer':
                if evalue <= hmmer_evalue:
                    fam, hmm_id = ref.split('_')[0], int(ref.split('_')[1])

                    fc_id = pc2fc_dict[fam][hmm_id]
                    func = fc2func_dict[fam][fc_id]
                    aln_coor = prot.split(':', 1)[1]
                    
                    # check if rdrp satisfy the requirements
                    try:
                        if func[0] == '*':
                            if evalue > rdrp_evalue or contig_len_dict[query] < rdrp_len:
                                continue
                    except KeyError:
                        pass

                    if query in group_hmmer_rst_dict[fam]:
                        group_hmmer_rst_dict[fam][query][prot] = fc_id
                    else:
                        group_hmmer_rst_dict[fam][query] = {prot: fc_id}

                    if func[0] == '*':
                        if evalue <= rdrp_evalue:
                            if contig_len_dict[query] >= rdrp_len:
                                fam_rdrp_set.add(fam)

    print(f"Identified families with RdRp:", fam_rdrp_set)

    # output the results step 2: calculate the completeness and conservation score for each family
    final_out = open(f"{temp_dir}/segvir.csv", 'w')
    final_out.write(f"family\tcontig\tlength\tmethod\tcoordinate\tevalue\ttaxon\tfunction\n")
    final_score_out = open(f"{temp_dir}/segvir.score.csv", 'w')
    final_score_out.write(f"family\tcompleteness\tcs\tref_cs\n")
    for fam in fam_rdrp_set:
        fc_set = set([])
        singleton_set = set([])
        temp_hmmer_aln_coor_dict = {}
        # print(fam, group_hmmer_rst_dict[fam], group_diamond_rst_dict[fam])
        for query, coor_dict in group_hmmer_rst_dict[fam].items():
            temp_hmmer_aln_coor_dict[query] = set([])
            for prot, fc in coor_dict.items():
                fc_set.add(fc)
                coor = prot.split(':', 1)[1]
                coor_start = min([int(coor.split(':')[0]), int(coor.split(':')[1])])
                coor_end = max([int(coor.split(':')[0]), int(coor.split(':')[1])])
                for i in range(coor_start, coor_end+1):
                    temp_hmmer_aln_coor_dict[query].add(i)

                try:
                    diamond_ref = coor2ref_dict[query][coor]
                    taxon = acc_name_dict[diamond_ref][1]
                except KeyError:
                    taxon = 'unknown'
                function = fc2func_dict[fam][fc][1]
                evalue = hmm_dict[query][prot][1]

                # discard the RdRp that doesn't meet the requirement
                if fc2func_dict[fam][fc][0] == '*':
                    if evalue > rdrp_evalue or contig_len_dict[query] < rdrp_len:
                        continue

                if fam in vir_single_seg_dict:
                    if taxon not in vir_single_seg_dict[fam]:
                        final_out.write(f"{fam}\t{query}\t{contig_len_dict[query]}\thmmer\t{coor}\t{evalue}\t{taxon}\t{function}\n")
                    else:
                        if evalue > 1e-20:
                            final_out.write(f"{fam}\t{query}\t{contig_len_dict[query]}\thmmer\t{coor}\t{evalue}\t{taxon}\t{function}\n")
                else:
                    final_out.write(f"{fam}\t{query}\t{contig_len_dict[query]}\thmmer\t{coor}\t{evalue}\t{taxon}\t{function}\n")

        for query, coor_dict in group_diamond_rst_dict[fam].items():
            for prot, ref in coor_dict.items():
                if ref in singleton2func_dict:
                    coor = prot.split(':', 1)[1]
                    coor_start = min([int(coor.split(':')[0]), int(coor.split(':')[1])])
                    coor_end = max([int(coor.split(':')[0]), int(coor.split(':')[1])])
                    diamond_coor_set = set([i for i in range(coor_start, coor_end+1)])
                    if query in temp_hmmer_aln_coor_dict:
                        if len(temp_hmmer_aln_coor_dict[query].intersection(diamond_coor_set))/ len(diamond_coor_set) <= 0.1:
                            singleton_set.add(ref)
                            taxon = acc_name_dict[ref][1]
                            function = singleton2func_dict[ref][1]
                            evalue = diamond_dict[query][prot][1]
                            final_out.write(f"{fam}\t{query}\t{contig_len_dict[query]}\tdiamond\t{coor}\t{evalue}\t{taxon}\t{function}\n")
                    else:
                        singleton_set.add(ref)
                        taxon = acc_name_dict[ref][1]
                        function = singleton2func_dict[ref][1]
                        evalue = diamond_dict[query][prot][1]
                        final_out.write(f"{fam}\t{query}\t{contig_len_dict[query]}\tdiamond\t{coor}\t{evalue}\t{taxon}\t{function}\n")

        print(fam, fc_set, singleton_set)

        if len(fc_set) == 0:
            final_score_out.write(f"{fam}\tnan\tnan\tnan\n")
            continue
                        
        # output the results step 3: family,completeness,conservation_score,contig,taxon,function
        # calcualte the completeness
        nearest_fc_set, ref_cs = set([]), 0.0
        temp_large_index = 0.0
        for ref_fc_set, score in fc_score_dict[fam].items():
            ref_fc_set = set([int(i) for i in ref_fc_set])
            set_inter = set.intersection(fc_set, ref_fc_set)
            set_union = set.union(fc_set, ref_fc_set)
            j_index = len(set_inter) * len(set_inter)/ len(fc_set)/ len(set_union)
            if j_index >= temp_large_index:
                nearest_fc_set = ref_fc_set
                ref_cs = score
                temp_large_index = j_index

        if len(nearest_fc_set) == 0:
            continue

        completeness = min(1.0, (len(fc_set) + len(singleton_set))/ len(nearest_fc_set))
        print(nearest_fc_set, ref_cs, completeness)

        # calculate the conservation score
        rdrp_score = []
        non_rdrp_score = []
        fc_set_weights = 0.0

        for rdrp_fc_id, fc_dict in fc2weight_dict[fam].items():
            temp_score = 0.0
            for fc in fc_set:
                if fc in fc_dict:
                    temp_score += fc_dict[fc]
            if rdrp_fc_id in fc_set:
                rdrp_score.append(temp_score)
            else:
                non_rdrp_score.append(temp_score + 1)
        
        if rdrp_score:
            fc_set_weights = max(rdrp_score)
        else:
            fc_set_weights = max(non_rdrp_score)
        print(fc_set_weights)

        final_score_out.write(f"{fam}\t{completeness}\t{fc_set_weights}\t{ref_cs}\n")

        # with open(f"{ref_dir}/segvir.csv", 'w') as final_out:
        #     final_out.write(f"family\tcontig\ttaxon\tfunction\n")
        
        # with open(f"{ref_dir}/segvir.score.csv", 'w') as final_out:
        #     final_out.write(f"family\tcompleteness\tconservation_score\n")

    final_out.close()
    final_score_out.close()
    segvir_df = pd.read_csv(f"{temp_dir}/segvir.csv", sep='\t')
    # segvir_df = segvir_df.sort_values('taxon').groupby('family').head()
    segvir_df = segvir_df.sort_values(['family'], ascending=True) \
    .groupby(['family'], sort=False) \
    .apply(lambda x: x.sort_values(['taxon'], ascending=True)) \
    .reset_index(drop=True)
    # segvir_df.sort_values(by=['family'], inplace=True)
    segvir_df.to_csv(f"{temp_dir}/segvir.csv", sep='\t', index=False)


def segvir_output(temp_dir, out_dir, outfmt=1):
    """Output the results of SegVir.
    """
    # load the parsed results
    contig_index = SeqIO.index(f"{temp_dir}/contig.clean.fna", 'fasta')
    
    host_coor_dict = {}
    host_df = pd.read_csv(f"{temp_dir}/candidate_blast.host.csv", sep='\t')
    for query, ref, aln_s, aln_e in zip(host_df['query'], host_df['reference'], host_df['aln_start'], host_df['aln_end']):
        if query not in host_coor_dict:
            host_coor_dict[query] = [(aln_s, aln_e)]
        else:
            host_coor_dict[query] = list(merge_range(host_coor_dict[query] + [(aln_s, aln_e)]))

    candidate_rst = f"{temp_dir}/candidate_result.csv"
    cd_df = pd.read_csv(candidate_rst, sep='\t')
    overlap_dict = {}
    for contig_id, overlap in zip(cd_df['query'], cd_df['host_ratio']):
        overlap_dict[contig_id] = overlap

    # output the contigs for each family
    shutil.copy(f"{temp_dir}/segvir.csv", f"{out_dir}")
    shutil.copy(f"{temp_dir}/segvir.score.csv", f"{out_dir}")
    segvir_df = pd.read_csv(f"{temp_dir}/segvir.csv", sep='\t')
    with open(f"{out_dir}/segvir.fna", 'w') as fna_out:
        for contig_id, family, coor in zip(segvir_df['contig'], segvir_df['family'], segvir_df['coordinate']):
            if contig_id in host_coor_dict and overlap_dict[contig_id] == 0:
                hs, he = min(host_coor_dict[contig_id][0]), max(host_coor_dict[contig_id][0])
                q_coor = [int(i) for i in coor.split(':')]
                qs, qe = min(q_coor), max(q_coor)
                contig_seq = contig_index[contig_id].seq
                if he < qs:
                    cut_seq = contig_seq[he: ]
                elif qe < hs:
                    cut_seq = contig_seq[: hs]
                fna_out.write(f">{contig_id} {family}\n{cut_seq}\n")
            else:
                fna_out.write(f">{contig_id} {family}\n{contig_index[contig_id].seq}\n")


if __name__ == "__main__":

    segvir_args = util.cmd_segvir()
    segvir_work_dir_path = os.getcwd()

    # input_contig = segvir_args.input
    # output_dir = segvir_args.outdir
    # ref_db_path = f"/home/tangxubo/Desktop/Seg_virus/SegVir_debug/test/ref_DB"
    # ref_db_dir = f"/home/tangxubo/Desktop/Seg_virus/SegVir_debug/test/hosts"
    # host_dir = f"/home/tangxubo/Desktop/Seg_virus/SegVir_debug/test/hosts"
    # temp_dir = f"/home/tangxubo/Desktop/Seg_virus/SegVir_debug/test/temp"
    # number_threads = 8

    # set input and output path
    input_contig = segvir_args.input
    output_dir = segvir_args.outdir
    # set the temp dir to save the temporary files
    if segvir_args.tempdir == 'temp':
        temp_dir = f"{output_dir}/temp"
    else:
        temp_dir = segvir_args.tempdir
    # load the path of reference database
    if segvir_args.database == 'segvir_db':
        ref_db_path = f"{segvir_work_dir_path}/segvir_db"
    else:
        ref_db_path = segvir_args.database
    # load the path of host genomes
    if segvir_args.host == None:
        host_path = None
    else:
        host_path = segvir_args.host
    # set the number of threads
    number_threads = segvir_args.thread

    for d in [output_dir, temp_dir]:
        if not os.path.exists(d):
            os.makedirs(d)
        # else:
        #     print(f"The output directory exists!")
        #     exit(0)

    if host_path:
        build_host_db(host_path=host_path, 
                    temp_dir=temp_dir)
        host_db_path=f"{temp_dir}/host"
    else:
        host_db_path=None
    
    run_segvir_aln(contig_path=input_contig, host_db_path=host_db_path, ref_dir=ref_db_path, temp_dir=temp_dir, 
                   threads=number_threads, len_thres=segvir_args.min_len)

    segvir_parse_rst(contig_path=input_contig, ref_dir=ref_db_path, temp_dir=temp_dir, 
                     blastp_evalue=segvir_args.blastp, hmmer_evalue=segvir_args.hmmer, host_ratio_thres=0.84,
                     rdrp_evalue=segvir_args.rdrp_evalue, rdrp_len=segvir_args.rdrp_len)

    segvir_output(temp_dir=temp_dir, 
                  out_dir=segvir_args.outdir)
