from Bio import SeqIO
import os
import pandas as pd
import util
import itertools
import collections


FAMILY_LIST = util.seg_family()


def filter_diamond(diamond_out_path, acc_taxon_dict, parse_out, bscore_thres=0.85):
    """Polish DIAMOND BLASTp output.
    """
    # Polish the diamond results, and parse out the detailed results
    # load the highest bit-score for each proteins
    prot_bs_dict = {}
    with open(diamond_out_path) as bo:
        for l in bo:
            records = l.split()
            prot_id = records[0]
            bscore = float(records[11])
            if prot_id in prot_bs_dict:
                if bscore >= prot_bs_dict[prot_id]:
                    prot_bs_dict[prot_id] = bscore
            else:
                prot_bs_dict[prot_id] = bscore
        
    # filter the results using the threshold of bit-score * 0.85
    with open(parse_out, 'w') as ps_out:
        ps_out.write(f"qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tnt_acc\tfamily\tref_name\n")
        with open(diamond_out_path) as bo:
            for l in bo:
                records = l.strip().split()
                prot_id = records[0]
                nt_id = records[0].split('_', 1)[1].split(':')[0]
                ref_id = records[1]
                bscore = float(records[11])
                if bscore >= prot_bs_dict[prot_id] * bscore_thres:
                    fam = acc_taxon_dict[ref_id][0]
                    virus = acc_taxon_dict[ref_id][1]
                    records.extend([nt_id, fam, virus])
                    records = '\t'.join(records)
                    ps_out.write(f"{records}\n")

    # sort the parsed results: grouped by nt_acc and sort by qseqid
    parse_out_df = pd.read_csv(parse_out, sep='\t')
    parse_out_df = parse_out_df.sort_values(['nt_acc'], ascending=True) \
    .groupby(['nt_acc'], sort=False) \
    .apply(lambda x: x.sort_values(['qseqid'], ascending=True)) \
    .reset_index(drop=True)
    parse_out_df.to_csv(parse_out, sep='\t', index=False)


def filter_hmmer(hmm_out_path, pc_fc_function_dict, parse_out, bscore_thres=0.85):
    """Parse HMMER output.
    """
    # parse the HMMER result, use the best hit, 
    # load the highest bit-score for each proteins
    prot_bs_dict = {}
    with open(hmm_out_path) as hop:
        for l in hop:
            if l[0] != '#' and l[0] != '':
                l = l.strip().split()
                prot_id = l[0]
                bscore = float(l[7])
                if prot_id in prot_bs_dict:
                    if bscore >= prot_bs_dict[prot_id]:
                        prot_bs_dict[prot_id] = bscore
                else:
                    prot_bs_dict[prot_id] = bscore
        
    # filter the results using the threshold of bit-score * 0.85
    with open(parse_out, 'w') as ps_out:
        ps_out.write(f"qseqid\thmm_id\tq_start\tq_end\tevalue\tbitscore\tnt_acc\tfamily\tFC\tfunction\n")
        with open(hmm_out_path) as ho:
            for l in ho:
                if l[0] != '#' and l[0] != '':
                    records = l.strip().split()
                    prot_id = records[0]
                    nt_id = records[0].split('_', 1)[1].split(':')[0]
                    pc_id = records[3]
                    bscore = float(records[7])
                    if bscore >= prot_bs_dict[prot_id] * bscore_thres:
                        fam = pc_id.split('_')[0]
                        fc_id, function = "*", '*'
                        if pc_id in pc_fc_function_dict:
                            fc_id = pc_fc_function_dict[pc_id][0]
                            function = pc_fc_function_dict[pc_id][1]
                        e_value = records[6]
                        q_start = records[17]
                        q_end = records[18]
                        records = '\t'.join([prot_id, pc_id, q_start, q_end, e_value, str(bscore), nt_id, fam, str(fc_id), function])
                        ps_out.write(f"{records}\n")

    # sort the parsed results: grouped by nt_acc and sort by qseqid
    parse_out_df = pd.read_csv(parse_out, sep='\t')
    parse_out_df = parse_out_df.sort_values(['nt_acc'], ascending=True) \
    .groupby(['nt_acc'], sort=False) \
    .apply(lambda x: x.sort_values(['qseqid'], ascending=True)) \
    .reset_index(drop=True)
    parse_out_df.to_csv(parse_out, sep='\t', index=False)


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


def parse_host_aln(contig_path, indent_seq_set, temp_dir):
    """Remove the host contigs, and identify contigs containing genes from hosts and virus.
    Due to the presence of nrEVEs and chimeric assembly contigs, it is necessary to identify and detect whether the detected virus may originate from nrEVE.
    diamond dictionary {query: [fam, segment, function, evalue, query start, query end, reference]}
    hmmer dictionary (query: [fam, segment, function, evalue, query start, query end, query_trans_id])
    """
    # build the index of contigs
    seq_index = SeqIO.index(contig_path, 'fasta')
    
    # save the identified contigs (DIAMOND & HMMER) in the temporary directory
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

    # analyze the blastn result (the best hit)
    # need to count all the alignment region for the host in order to distinguish the chimeric assembly contigs and nrEVE
    contig_aln_dict = {}
    out_aln_coor_dict = {}
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
                out_aln_coor_dict[query] = [(qstart, qend)]
            else:
                if reference == contig_aln_dict[query][5]:
                    # if ident >= 80 and evalue <= 1e-5:
                    contig_aln_dict[query][6] = list(merge_range(contig_aln_dict[query][6] + [(qstart, qend)]))
                    out_aln_coor_dict[query] = list(merge_range(out_aln_coor_dict[query] + [(qstart, qend)]))

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
    
    return out_aln_coor_dict


def tuple_list2set(tuple_list):
    """Convert the list of tuples to a set of numbers.
    """
    out_set = set([])
    for (t1, t2) in tuple_list:
        for i in range(t1, t2+1):
            out_set.add(i)
    return out_set


def prot2taxon(prot_df, mode='best'):
    """Determine the alignemnt based on the bitsocre.
    Return the results based on the bitsocre..
    Input: {protein: [[fam, biscore]]}
    'best' output: {protein: fam}
    'all' output: {protein: [fam1, fam2, ...]}
    """
    prot_aln_dict = {}
    for qseqid, fam, bscore in zip(prot_df['qseqid'], prot_df['family'], prot_df['bitscore']):
        if qseqid in prot_aln_dict:
            prot_aln_dict[qseqid].append([fam, bscore])
        else:
            prot_aln_dict[qseqid] = [[fam, bscore]]

    prot2fam_score_dict = {}
    for prot, fs_list in prot_aln_dict.items():
        for (fam, score) in fs_list:
            if prot in prot2fam_score_dict:
                if fam in prot2fam_score_dict[prot]:
                    prot2fam_score_dict[prot][fam] += score
                else:
                    prot2fam_score_dict[prot][fam] = score
            else:
                prot2fam_score_dict[prot] = {}
                prot2fam_score_dict[prot][fam] = score
    
    prot2taxon_dict = {}
    if mode == 'best':
        for prot, fam_score_dict in prot2fam_score_dict.items():
            fam_score_dict = dict(sorted(fam_score_dict.items(), key=lambda item: item[1], reverse=True))
            prot2taxon_dict[prot] = list(fam_score_dict.keys())[0]
    elif mode == 'all':
        for prot, fam_score_dict in prot2fam_score_dict.items():
            prot2taxon_dict[prot] = list(fam_score_dict.keys())
    else:
        print('Unknown mode.')
        exit(0)

    return prot2taxon_dict


def segvir_parse_rst(ref_dir, temp_dir, blastp_evalue=1e-5, hmmer_evalue=1e-5, host_ratio_thres=0.84,
                     rdrp_evalue=1e-10, rdrp_len=900, vote_thres=0.85, len_thres=300):
    """Parse the alignment and output the final result.
    """
    contig_path=f"{temp_dir}/contig.clean.fna"
    
    # load the index of all contigs
    contig_seq_index = {}
    for s in SeqIO.parse(contig_path, 'fasta'):
        contig_seq_index[s.id] = s
    contig_len_dict = {s: len(contig_seq_index[s]) for s in contig_seq_index}

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

    # load the dictionary of pc: [fc, funtion]
    pc2fc_func_dict = {}
    fc_info_df = pd.read_csv(f"{ref_dir}/pc2fc.tsv", sep='\t')
    for family, pc, fc, anno_s, func in zip(fc_info_df['family'], fc_info_df['pc'], fc_info_df['fc'], fc_info_df['anno_state'], fc_info_df['function']):
        if pd.isna(func):
            func = 'unknown'
        pc2fc_func_dict[f'{family}_{pc}'] = [fc, func]


    # parse the results of Diamond and HMMER
    filter_diamond(diamond_out_path=f"{temp_dir}/contig.diamond", 
                   acc_taxon_dict=acc_name_dict, 
                   parse_out=f"{temp_dir}/contig.diamond.filter", 
                   bscore_thres=vote_thres)
    
    filter_hmmer(hmm_out_path=f"{temp_dir}/hmm.tblout", 
                 pc_fc_function_dict=pc2fc_func_dict, 
                 parse_out=f"{temp_dir}/hmm.tblout.filter", 
                 bscore_thres=vote_thres)

    # remove contimanition from rna virus using HMM out (the contig will be removed only if it contigs FP protein)
    hmm_rna_conti_set = set()
    hmm_filter_df = pd.read_csv(f"{temp_dir}/hmm.tblout.filter", sep='\t')
    # temp_dict = {}
    # for qseqid, fam, bscore in zip(hmm_filter_df['qseqid'], hmm_filter_df['family'], hmm_filter_df['bitscore']):
    #     if qseqid in temp_dict:
    #         temp_dict[qseqid].append([fam, bscore])
    #     else:
    #         temp_dict[qseqid] = [[fam, bscore]]
    # prot2taxon_dict = prot2taxon(prot_aln_dict=temp_dict, mode='best')
    prot2taxon_dict = prot2taxon(prot_df=hmm_filter_df, mode='best')
    for prot_id, taxon in prot2taxon_dict.items():
        if taxon == "FP":
            hmm_rna_conti_set.add(prot_id.split('_', 1)[1].split(':')[0])
    # print(f"The contigs from non-segmented RNA virues: {hmm_rna_conti_set}")

    # remove ambiguous contig using diamond out when contig encoded multiple proteins
    # the proteins of the contig should be mainly from the same family (> 50%)
    diamond_amb_set = set()
    diamond_filter_df = pd.read_csv(f"{temp_dir}/contig.diamond.filter", sep='\t')
    # temp_dict = {}
    # for qseqid, fam, bscore in zip(diamond_filter_df['qseqid'], diamond_filter_df['family'], diamond_filter_df['bitscore']):
    #     if qseqid in temp_dict:
    #         temp_dict[qseqid].append([fam, bscore])
    #     else:
    #         temp_dict[qseqid] = [[fam, bscore]]
    # prot2taxon_dict = prot2taxon(prot_aln_dict=temp_dict, mode='all')
    prot2taxon_dict = prot2taxon(prot_df=diamond_filter_df, mode='all')
    nt2taxon_dict = {}
    for prot_id, taxon_list in prot2taxon_dict.items():
        nt_id = prot_id.split('_', 1)[1].split(':')[0]
        if nt_id in nt2taxon_dict:
            nt2taxon_dict[nt_id].append(taxon_list)
    for nt, prot_taxon_list in nt2taxon_dict.items():
        all_taxon_list = list(itertools.chain.from_iterable(prot_taxon_list))
        counter = collections.Counter(all_taxon_list)
        [(taxon, freq)] = counter.counter.most_common(1)
        if freq/ len(prot_taxon_list) <= 0.5:
            diamond_amb_set.add(nt)
    # print(f"The contigs contains ambigious taxons: {diamond_amb_set}")

    # generate the results without rna virus and cmbigiuos contigs
    temp_set = set.union(hmm_rna_conti_set, diamond_amb_set)
    hmm_filter_df.drop(hmm_filter_df[hmm_filter_df['nt_acc'].isin(temp_set)].index, inplace=True)
    diamond_filter_df.drop(diamond_filter_df[diamond_filter_df['nt_acc'].isin(temp_set)].index, inplace=True)
    hmm_filter_df.to_csv(f"{temp_dir}/hmm.tblout.filter.no_rna.no_amb", sep='\t', index=False)
    diamond_filter_df.to_csv(f"{temp_dir}/contig.diamond.filter.no_rna.no_amb", sep='\t', index=False)

    # load the alignment regions of segmented RNA virus and hosts and remove the contigs if overlap >= 85%
    # diamond_ident_set = set(pd.read_csv(f"{temp_dir}/contig.diamond.filter", sep='\t')['qseqid'].values)
    # hmmer_ident_set = set(pd.read_csv(f"{temp_dir}/hmm.tblout.filter", sep='\t')['qseqid'].values)
    aln_host_coor_dict = parse_host_aln(contig_path=contig_path, 
                                        indent_seq_set=set.union(set(diamond_filter_df['nt_acc'].values), set(hmm_filter_df['nt_acc'].values)), 
                                        temp_dir=temp_dir)
    aln_vir_coor_dict = {}
    for query, prot_id in zip(diamond_filter_df['nt_acc'], diamond_filter_df['qseqid']):
        orf_s = min(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
        orf_e = max(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
        if query not in aln_vir_coor_dict:
            aln_vir_coor_dict[query] = [(orf_s, orf_e)]
        else:
            aln_vir_coor_dict[query] = list(merge_range(aln_vir_coor_dict[query] + [(orf_s, orf_e)]))
    for query, prot_id in zip(hmm_filter_df['nt_acc'], hmm_filter_df['qseqid']):
        orf_s = min(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
        orf_e = max(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
        if query not in aln_vir_coor_dict:
            aln_vir_coor_dict[query] = [(orf_s, orf_e)]
        else:
            aln_vir_coor_dict[query] = list(merge_range(aln_vir_coor_dict[query] + [(orf_s, orf_e)]))

    with open(f"{temp_dir}/vir.aln_coor.csv", 'w') as coor_out:
        for query, coor_list in aln_vir_coor_dict.items():
            for coor in coor_list:
                coor_out.write(f"{query}\t{coor[0]}:{coor[1]}\n")

    host_cont_set = set([])
    with open(f"{temp_dir}/host_overlap.ratio.csv", 'w') as overlap_out:
        for nt_acc, aln_coor_list in aln_host_coor_dict.items():
            host_coor_set = tuple_list2set(aln_coor_list)
            vir_coor_set = tuple_list2set(aln_vir_coor_dict[nt_acc])
            overlap_ratio = len(host_coor_set.intersection(vir_coor_set))/ len(vir_coor_set)
            overlap_out.write(f"{nt_acc}\t{overlap_ratio}\n")
            if overlap_ratio >= host_ratio_thres:
                host_cont_set.add(nt_acc)
    # print(f"The ambigious contigs contains host genes: {host_cont_set}")
    hmm_filter_df.drop(hmm_filter_df[hmm_filter_df['nt_acc'].isin(host_cont_set)].index, inplace=True)
    diamond_filter_df.drop(diamond_filter_df[diamond_filter_df['nt_acc'].isin(host_cont_set)].index, inplace=True)
    hmm_filter_df.to_csv(f"{temp_dir}/hmm.tblout.filter.no_rna.no_amb.no_host", sep='\t', index=False)
    diamond_filter_df.to_csv(f"{temp_dir}/contig.diamond.filter.no_rna.no_amb.no_host", sep='\t', index=False)

    # integrate the results of diamond and HMMER
    # load the singleton protein information (protein: function)
    singleton2func_dict = {}    # singleton: {anno_s, func}
    singleton_df = pd.read_csv(f"{ref_dir}/prot2func.csv", sep='\t')
    for prot_id, func in zip(singleton_df['protein_id'], singleton_df['function']):
        anno_s = '-'
        # if any(i in func.lower() for i in ['rdrp', 'polymerase']):
        #     anno_s = '*'
        singleton2func_dict[prot_id] = [anno_s, func]
    # load the set of single-segment viruses
    vir_single_set = set([])
    single_vir_df = pd.read_csv(f"{ref_dir}/virus.single_seg.csv", sep='\t')
    for vir_name in single_vir_df['virus']:
        vir_single_set.add(vir_name.lower())
    
    # Strategy 1 (best taxon based on KNN)
    hmm_best_dict = {}  # prot_id: bitscore, evalue, fam, fc, function
    prot2taxon_dict = prot2taxon(prot_df=hmm_filter_df, mode='best')
    for qseqid, fam, bscore, evalue, fc_id, function in zip(hmm_filter_df['qseqid'], hmm_filter_df['family'], hmm_filter_df['bitscore'], hmm_filter_df['evalue'], hmm_filter_df['FC'], hmm_filter_df['function']):
        if qseqid in prot2taxon_dict:
            if prot2taxon_dict[qseqid] == fam:
                if qseqid in hmm_best_dict:
                    if bscore >= hmm_best_dict[qseqid][0]:
                        hmm_best_dict[qseqid] = [bscore, evalue, fam, fc_id, function]
                else:
                    hmm_best_dict[qseqid] = [bscore, evalue, fam, fc_id, function]

    diamond_best_dict = {}  # prot_id: bitscore, evalue, fam, virus_name, ref_id
    prot2taxon_dict = prot2taxon(prot_df=diamond_filter_df, mode='best')
    for qseqid, ref_id, fam, bscore, evalue, virus_name in zip(diamond_filter_df['qseqid'], diamond_filter_df['sseqid'], diamond_filter_df['family'], diamond_filter_df['bitscore'], diamond_filter_df['evalue'], diamond_filter_df['ref_name']):
        if qseqid in prot2taxon_dict:
            if prot2taxon_dict[qseqid] == fam:
                if qseqid in diamond_best_dict:
                    if bscore >= diamond_best_dict[qseqid][0]:
                        diamond_best_dict[qseqid] = [bscore, evalue, fam, virus_name, ref_id]
                else:
                    diamond_best_dict[qseqid] = [bscore, evalue, fam, virus_name, ref_id]

    with open(f"{temp_dir}/merge_diamond_hmmer.csv", 'w') as rst_out:
        rst_out.write(f"family\tcontig\tlength\tgene_coor\tblastp_e\thmmer_e\tref_name\tfunction\tfc\n")
        diamond_set = set(list(diamond_best_dict.keys()))
        hmm_set = set(list(hmm_best_dict.keys()))
        for prot_id in set.intersection(diamond_set, hmm_set):
            if hmm_best_dict[prot_id][2] == diamond_best_dict[prot_id][2]:
                fam = diamond_best_dict[prot_id][2]
                contig_id = prot_id.split('_', 1)[1].split(':')[0]
                contig_len = contig_len_dict[contig_id]
                orf_s = min(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
                orf_e = max(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
                blastp_e = diamond_best_dict[prot_id][1]
                hmmer_e = hmm_best_dict[prot_id][1]
                # filter the result by e-values of BLASTp and HMMER
                if blastp_e > blastp_evalue and hmmer_e > hmmer_evalue:
                    continue
                virus_name = diamond_best_dict[prot_id][3]
                func = hmm_best_dict[prot_id][4]
                # filter the RdRps
                if any(i in func.lower() for i in ['rdrp', 'polymerase']):
                    if blastp_e > rdrp_evalue and hmmer_e > rdrp_evalue:
                        continue
                    else:
                        if contig_len < rdrp_len:
                            continue
                
                fc_id = hmm_best_dict[prot_id][3]
                if virus_name.lower() in vir_single_set:
                    if blastp_e <= 1e-25:
                        continue
                rst_out.write(f"{fam}\t{contig_id}\t{contig_len}\t{orf_s}:{orf_e}\t{blastp_e}\t{hmmer_e}\t{virus_name}\t{func}\t{fc_id}\n")
        
        for prot_id in diamond_set - hmm_set:
            fam = diamond_best_dict[prot_id][2]
            contig_id = prot_id.split('_', 1)[1].split(':')[0]
            contig_len = contig_len_dict[contig_id]
            orf_s = min(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
            orf_e = max(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
            blastp_e = diamond_best_dict[prot_id][1]
            if blastp_e > blastp_evalue:
                continue
            hmmer_e = '*'
            virus_name = diamond_best_dict[prot_id][3]
            func = singleton2func_dict[diamond_best_dict[prot_id][4]][1]
            # filter the RdRps
            if any(i in func.lower() for i in ['rdrp', 'polymerase']):
                if blastp_e > rdrp_evalue:
                    continue
                else:
                    if contig_len < rdrp_len:
                        continue
            fc_id = '*'
            if virus_name.lower() in vir_single_set:
                if blastp_e <= 1e-25:
                    continue
            rst_out.write(f"{fam}\t{contig_id}\t{contig_len}\t{orf_s}:{orf_e}\t{blastp_e}\t{hmmer_e}\t{virus_name}\t{func}\t{fc_id}\n")

        for prot_id in hmm_set - diamond_set:
            fam = hmm_best_dict[prot_id][2]
            contig_id = prot_id.split('_', 1)[1].split(':')[0]
            contig_len = contig_len_dict[contig_id]
            orf_s = min(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
            orf_e = max(int(prot_id.split(':')[-2]), int(prot_id.split(':')[-1]))
            blastp_e = '*'
            hmmer_e = hmm_best_dict[prot_id][1]
            if hmmer_e > hmmer_evalue:
                continue
            virus_name = '*'
            func = hmm_best_dict[prot_id][4]
            # filter the RdRps
            if any(i in func.lower() for i in ['rdrp', 'polymerase']):
                if hmmer_e > rdrp_evalue:
                    continue
                else:
                    if contig_len < rdrp_len:
                        continue
            fc_id = hmm_best_dict[prot_id][3]
            rst_out.write(f"{fam}\t{contig_id}\t{contig_len}\t{orf_s}:{orf_e}\t{blastp_e}\t{hmmer_e}\t{virus_name}\t{func}\t{fc_id}\n")
    
    # filter out the families without RdRps
    rdrp_fam_set = set()
    merge_df = pd.read_csv(f"{temp_dir}/merge_diamond_hmmer.csv", sep='\t')
    for fam, func in zip(merge_df['family'], merge_df['function']):
        if any(i in func.lower() for i in ['rdrp', 'polymerase']):
            rdrp_fam_set.add(fam)
    merge_df = merge_df[merge_df['family'].isin(rdrp_fam_set)]
    merge_df.to_csv(f"{temp_dir}/merge_diamond_hmmer.rdrp.csv", sep='\t', index=False)

    # filter out the contigs with shorter length
    merge_df = merge_df[merge_df['length'] >= len_thres]
    
    # output the detailed information of SegVir's outputs
    merge_df = merge_df.sort_values(['family'], ascending=True) \
    .groupby(['family'], sort=False) \
    .apply(lambda x: x.sort_values(['ref_name'], ascending=True)) \
    .reset_index(drop=True)
    merge_df.to_csv(f"{output_dir}/segvir.csv", sep='\t', index=False)

    # calcualte the estimated completeness and conservation scores
    # load the fc weights
    fc2weight_dict = {f: {} for f in FAMILY_LIST}  # family: {rdrp_fc_id: {fc_id: weights}}
    weight_df = pd.read_csv(f"{ref_dir}/fc_weight.csv", sep='\t')
    for fam, rdrp_fc, fc, weight in zip(weight_df['family'], weight_df['rdrp_fc'], weight_df['fc'], weight_df['weight']):
        if rdrp_fc in fc2weight_dict[fam]:
            fc2weight_dict[fam][rdrp_fc][fc] = weight
        else:
            fc2weight_dict[fam][rdrp_fc] = {fc: weight}
    # load the reference fc sets and conservation score
    fc_score_dict = {f: {} for f in FAMILY_LIST}
    fc_score_df = pd.read_csv(f"{ref_dir}/ref_fc_set_cs.tsv", sep='\t')
    for fam, fc_set, score in zip(fc_score_df['family'], fc_score_df['fc_set'], fc_score_df['score']):
        fc_set = frozenset([i for i in fc_set.split(',')])
        if fc_set not in fc_score_dict[fam]:
            fc_score_dict[fam][fc_set] = score

    fam_fc_set_dict = {}
    for fam, fc in zip(merge_df['family'], merge_df['fc']):
        if fc != '*':
            if fam in fam_fc_set_dict:
                fam_fc_set_dict[fam].add(fc)
            else:
                fam_fc_set_dict[fam] = set([fc])

    fam_singleton_dict = {}
    for fam, contig, fc in zip(merge_df['family'], merge_df['contig'], merge_df['fc']):
        if fc == '*':
            if fam in fam_singleton_dict:
                fam_singleton_dict[fam].add(contig)
            else:
                fam_singleton_dict[fam] = set([contig])

    with open(f"{output_dir}/segvir.score.csv", 'w') as final_score_out:
        for fam, fc_set in fam_fc_set_dict.items():
            if len(fc_set) == 0:
                final_score_out.write(f"{fam}\tnan\tnan\tnan\n")
                continue

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
            
            num_singleton = 0
            if fam in fam_singleton_dict:
                num_singleton = len(fam_singleton_dict[fam])
            completeness = min(1.0, (len(fc_set) + num_singleton)/ len(nearest_fc_set))
            # print(nearest_fc_set, ref_cs, completeness)

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
            # print(fc_set_weights)

            final_score_out.write(f"{fam}\t{completeness}\t{fc_set_weights}\t{ref_cs}\n")


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

    candidate_rst = f"{temp_dir}/host_overlap.ratio.csv"
    overlap_df = pd.read_csv(candidate_rst, sep='\t', header=None)
    overlap_dict = {}
    for contig_id, overlap in zip(overlap_df[0], overlap_df[1]):
        overlap_dict[contig_id] = overlap

    # load the results of SegVir
    segvir_df = pd.read_csv(f"{out_dir}/segvir.csv", sep='\t')

    # remove the chimeric regions from contigs
    seq_region_dict = {}
    for contig_id, coor in zip(segvir_df['contig'], segvir_df['gene_coor']):
        q_coor = [int(i) for i in coor.split(':')]
        if contig_id in seq_region_dict:
            seq_region_dict[contig_id].extend(q_coor)
        else:
            seq_region_dict[contig_id] = q_coor

    seq_dict = {}
    for contig_id, family, coor, function in zip(segvir_df['contig'], segvir_df['family'], segvir_df['gene_coor'], segvir_df['function']):
        if contig_id not in seq_dict:
            if contig_id in host_coor_dict and overlap_dict[contig_id] == 0:
                hs, he = min(host_coor_dict[contig_id][0]), max(host_coor_dict[contig_id][0])
                qs, qe = min(seq_region_dict[contig_id]), max(seq_region_dict[contig_id])
                contig_seq = contig_index[contig_id].seq
                print(contig_id, hs, he, qs, qe)
                if he <= qs:
                    cut_seq = contig_seq[he: ]
                elif qe < hs:
                    cut_seq = contig_seq[: hs]
                seq_dict[contig_id] = cut_seq
            else:
                seq_dict[contig_id] = str(contig_index[contig_id].seq)

    # save the results in dictionary
    fam_dict = {}
    for contig_id, family, coor, function in zip(segvir_df['contig'], segvir_df['family'], segvir_df['gene_coor'], segvir_df['function']):
        anno_s = '-'
        if any(i in function.lower() for i in ['rdrp', 'polymerase']):
            anno_s = 'rdrp'
        else:
            anno_s = 'nonrdrp'
        
        if family in fam_dict:
            fam_dict[family][anno_s].add(contig_id)
        else:
            fam_dict[family] = {'rdrp': set([]), 'nonrdrp': set([])}
            fam_dict[family][anno_s] = set([contig_id])

    # remove the duplicated sequence in non-RdRp set
    for fam in fam_dict:
        fam_dict[fam]['nonrdrp'] = fam_dict[fam]['nonrdrp'] - fam_dict[fam]['rdrp']

    if outfmt == 1:
        with open(f"{out_dir}/segvir.fna", 'w') as fna_out:
            for fam, rdrp_dict in fam_dict.items():
                for sid in rdrp_dict['rdrp']:
                    fna_out.write(f">{sid} {fam}\n{seq_dict[sid]}\n")
                for sid in rdrp_dict['nonrdrp']:
                    fna_out.write(f">{sid} {fam}\n{seq_dict[sid]}\n")

    elif outfmt == 2:
        if not os.path.exists(f"{out_dir}/SegVir_seq"):
            os.mkdir(f"{out_dir}/SegVir_seq")
        for fam, rdrp_dict in fam_dict.items():
            with open(f"{out_dir}/SegVir_seq/{fam}.fna", 'w') as fna_out:
                for sid in rdrp_dict['rdrp']:
                    fna_out.write(f">{sid}\n{seq_dict[sid]}\n")
                for sid in rdrp_dict['nonrdrp']:
                    fna_out.write(f">{sid}\n{seq_dict[sid]}\n")
    
    elif outfmt == 3:
        if not os.path.exists(f"{out_dir}/SegVir_seq"):
            os.mkdir(f"{out_dir}/SegVir_seq")
        for fam, rdrp_dict in fam_dict.items():
            if len(rdrp_dict['rdrp']) > 0:
                with open(f"{out_dir}/SegVir_seq/{fam}.rdrp.fna", 'w') as fna_out:
                    for sid in rdrp_dict['rdrp']:
                        fna_out.write(f">{sid}\n{seq_dict[sid]}\n")
            
            if len(rdrp_dict['nonrdrp']) > 0:
                with open(f"{out_dir}/SegVir_seq/{fam}.nonrdrp.fna", 'w') as fna_out:
                    for sid in rdrp_dict['nonrdrp']:
                        fna_out.write(f">{sid}\n{seq_dict[sid]}\n")


if __name__ == "__main__":

    segvir_args = util.cmd_segvir_polish()
    segvir_work_dir_path = os.getcwd()

    # load the path of reference database
    if segvir_args.database == 'segvir_db':
        ref_db_path = f"{segvir_work_dir_path}/segvir_db"
    else:
        ref_db_path = segvir_args.database

    # check if the output directory exists
    output_dir = segvir_args.outdir
    overwrite_bool = segvir_args.overwrite
    if os.path.exists(output_dir):
        if os.path.exists(f"{output_dir}/segvir.csv"):
            if not overwrite_bool:
                print(f"The output results exists. (Please set `--overwrite` as True to overwrite the results.)")
                exit(0)
    else:
        os.makedirs(output_dir)
    
    segvir_parse_rst(ref_dir=ref_db_path, temp_dir=segvir_args.tempdir, 
                     blastp_evalue=segvir_args.blastp, hmmer_evalue=segvir_args.hmmer, 
                     host_ratio_thres=0.5, rdrp_evalue=segvir_args.rdrp_evalue, 
                     rdrp_len=segvir_args.rdrp_len, vote_thres=segvir_args.vote_thres,
                     len_thres=segvir_args.min_len)

    segvir_output(temp_dir=segvir_args.tempdir, out_dir=segvir_args.outdir,
                  outfmt=segvir_args.outfmt)
