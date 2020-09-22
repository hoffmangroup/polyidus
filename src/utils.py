from datetime import datetime
import numpy as np
import os
import pandas as pd
import pysam
import subprocess


CIGAR_DEFS = {0: "Match", 1: "Insertion", 2: "Deletion",
              3: "RefSkip", 4: "SoftClip", 5: "HardClip",
              6: "Pad", 7: "Equal", 8: "Diff",
              9: "Back", 10: "Unmapped"}


def get_chimer(chimerpath, minreads, minscore):
    chimerdf = pd.read_csv(chimerpath, sep="\t")
    print(chimerdf.head())
    chimerdf = chimerdf[chimerdf["NumberReads"] > minreads]
    chimerdf = chimerdf[chimerdf["Normalized.Score"] > minscore]
    return chimerdf


def get_virus_dict(viralbam, readnames):
    bam = pysam.AlignmentFile(viralbam, "rb")
    dict_viral = {}
    for read in bam:
        if read.qname in readnames or read.qname + ".1" in readnames:
            if len(read.cigar) > 1 and not read.is_secondary:
                readname = read.qname + ".1"
                if read.is_read2:
                    readname = read.qname + ".2"
                if readname in dict_viral.keys():
                    print("Error!")
                dict_viral[readname] = read
    return dict_viral


def get_virus_info(dict_viral, readname, window):
    read1 = dict_viral.get(readname + ".1", "")
    read2 = dict_viral.get(readname + ".2", "")
    try:
        chromlen = read1.header.lengths[0]
    except Exception:
        chromlen = 7904
    reads = []
    if not isinstance(read1, str):
        reads.append(read1)
    if not isinstance(read2, str):
        reads.append(read2)
    sts = []
    ends = []
    for read in reads:
        read_positions = read.get_reference_positions(full_length=True)
        np_none = np.where(
            [np.logical_not(np.isscalar(each)) for
             each in read_positions])[0]
        if len(read.cigar) == 2:
            if 0 in np_none:
                adst = read_positions[max(np_none) + 1]
                adend = adst + window
            else:
                adend = read_positions[min(np_none) - 1]
                adst = adend - window
            if adend not in ends and adst not in sts:
                ends.append(adend)
                sts.append(adst)
    if len(sts) == 0:
        sts = [0]
        ends = [0]
        chromvir = "NA"
    else:
        chromvir = read.reference_name
    return sts[0], ends[0], chromvir, chromlen


def get_host_pos(viralbam, hostbam, chromhost, start_host,
                 readnames, window=10000):
    dict_viral = get_virus_dict(viralbam, readnames)
    bam = pysam.AlignmentFile(hostbam, "rb")
    dict_supports = {}
    st_hosts = []
    end_hosts = []
    host_integs = []
    st_virs = []
    end_virs = []
    viral_integs = []
    vir_chroms = []
    new_read_names = []
    dict_readnames = {}
    # Suggested by johnstonmj to make sure not exceeding chrom bounds
    chromsize = bam.header.get_reference_length(chromhost)
    lower_bound = max([0, start_host - 2000])
    upper_bound = min([start_host + 2000, chromsize])
    # Iterate through the reads within the bounds
    for read in bam.fetch(chromhost, lower_bound, upper_bound):
        if read.qname in readnames or read.qname + ".1" in readnames:
            if len(read.cigar) > 1:
                read_positions = read.positions
                two_cigars = len(read.cigar) == 2
                # two_cigars = len(read.cigar) >= 2
                if len(read_positions) < len(read.seq) and two_cigars:
                    virst, virend, virchrom, chromlen = get_virus_info(
                        dict_viral, read.qname, window)
                    if virst != virend:
                        new_read_names.append(read.qname)
                        read_positions = read.get_reference_positions(
                            full_length=True)
                        np_none = np.where(
                            [np.logical_not(np.isscalar(each)) for
                             each in read_positions])[0]
                        len_poses = len(read_positions)
                        if 0 in np_none and (len_poses - 1) not in np_none:
                            host_integ = read_positions[max(np_none) + 1]
                            adst = read_positions[max(np_none) + 1]
                            adend = adst + window
                        elif 0 in np_none:
                            host_integ = read_positions[np_none[1] - 1]
                            adend = read_positions[np_none[1] - 1]
                            adst = adend - window
                        else:
                            host_integ = read_positions[min(np_none) - 1]
                            adend = read_positions[min(np_none) - 1]
                            adst = adend - window
                        dict_supports[host_integ] = \
                            dict_supports.get(host_integ, 0) + 1
                        # Keep track of read names used for each integration
                        list_curreadnames = dict_readnames.get(host_integ, [])
                        list_curreadnames.append(read.qname)
                        dict_readnames[host_integ] = list_curreadnames
                        if host_integ not in host_integs:
                            host_integs.append(host_integ)
                            end_hosts.append(adend)
                            viral_integ = virst
                            if virst < 0 or virst > 7904:
                                viral_integ = virend
                            viral_integs.append(viral_integ)
                            st_hosts.append(adst)
                            st_virs.append(virst)
                            end_virs.append(virend)
                            vir_chroms.append(virchrom)
    if len(vir_chroms) == 0:
        vir_chroms = ["NA"]
        chromlen = 7904
    dict_out = {"ChromHost": [chromhost for each in range(len(st_hosts))],
                "Starts": st_hosts,
                "Ends": end_hosts,
                "VirusStarts": st_virs,
                "VirusEnds": end_virs,
                "ChromVirus": [vir_chroms[0] for each in range(len(st_hosts))],
                "VirusChromLen": chromlen,
                "IntegrationHosts": host_integs,
                "IntegrationVirus": viral_integs,
                "NumberReads": [dict_supports.get(each, 0)
                                for each in host_integs],
                "ReadNameDict": dict_readnames}
    return dict_out


def get_virus_seq(fastqvirusobj, chrom, st, end, chromlen):
    sts = []
    ends = []
    seq = ""
    if st < 0:
        virusst = 0
        virusend = end
        sts.append(0)
        ends.append(end)
        adseq = fastqvirusobj.fetch(chrom, 0, end)
        seq = seq + adseq
        sts.append(chromlen)
        ends.append(end)
        adseq = fastqvirusobj.fetch(chrom, end, chromlen)
        seq = adseq + seq
    if end > chromlen:
        virusst = st
        virusend = chromlen
        sts.append(st)
        ends.append(chromlen)
        adseq = fastqvirusobj.fetch(chrom, st, chromlen)
        seq = seq + adseq
        sts.append(0)
        ends.append(st)
        adseq = fastqvirusobj.fetch(chrom, 0, st)
        seq = seq + adseq
    if end < 0:
        virusst = st
        virusend = 0
        sts.append(st)
        ends.append(0)
        adseq = fastqvirusobj.fetch(chrom, 0, st)
        seq = seq + adseq[::-1]
        sts.append(chromlen)
        ends.append(st)
        adseq = fastqvirusobj.fetch(chrom, st, chromlen)
        seq = seq + adseq[::-1]
    return seq, virusst, virusend


def make_rev_comp(seq):
    seq_dict = {"A": "T",
                "T": "A",
                "G": "C",
                "C": "G",
                "N": "N"}
    revcomp = ["A" for each in seq]
    i = 0
    for each_nuc in seq[::-1]:
        try:
            revcomp[i] = seq_dict[each_nuc.decode()]
        except:
            revcomp[i] = seq_dict[each_nuc]
        i = i + 1
    rev_str = "".join(revcomp)
    return revcomp, rev_str


def get_read_dict(bam_path):
    sam_link = pysam.AlignmentFile(bam_path, "rb")
    # order_chroms = []
    dict_reads = {}
    for read in sam_link:
        if not read.is_secondary:
            read_name = read.qname
            dict_reads = parse_read(read, dict_reads, read_name)
    return dict_reads


def parse_read(read, dict_reads, read_name):
    keep_read = False
    has_softclip = False
    cigar_stat = {}
    if read.cigartuples is not None:
        idx_cigar = 0
        list_names = []
        length_match = 0
        for cigar_tuple in read.cigartuples:
            ad_name = CIGAR_DEFS[cigar_tuple[0]]
            if ad_name in list_names:
                ad_name = ad_name + ".2"
            cigar_stat[idx_cigar] = {
                "Name": ad_name,
                "Length": cigar_tuple[1]}
            if ad_name == "Match":
                length_match = cigar_tuple[1]
            idx_cigar += 1
            list_names.append(ad_name)
        if "SoftClip.2" not in list_names:
            keep_read = True
            seq_pos = read.get_reference_positions()
            start, end = [min(seq_pos), max(seq_pos)]
            mapping_frac = length_match / (end - start + 1)
            if "S" in read.cigarstring:
                has_softclip = True
            if read_name in dict_reads.keys():
                prev_frac = dict_reads[read_name]["Mapping fraction"]
                if mapping_frac < prev_frac:
                    keep_read = False
            if keep_read:
                dict_reads[read_name] = {
                    "CIGAR": cigar_stat,
                    "CIGARSTR": read.cigarstring,
                    "Chrom": read.reference_name,
                    "Sequence": read.seq,
                    "Start": start,
                    "IsPrimary": np.logical_not(read.is_secondary),
                    "Negative strand": read.is_reverse,
                    "End": end,
                    "Mapping fraction": mapping_frac,
                    "HasSoftClip": has_softclip,
                    "ReadName": read_name}
    return dict_reads


def get_paired_end_cigar(read1, read2):
    CIGAR_DEFS = {0: "Match", 1: "Insertion", 2: "Deletion",
                  3: "RefSkip", 4: "SoftClip", 5: "HardClip",
                  6: "Pad", 7: "Equal", 8: "Diff",
                  9: "Back", 10: "Unmapped"}
    cigar_stat = {}
    idx_cigar = 0
    list_names = []
    if read1.cigarstring is None:
        cigar_stat[idx_cigar] = {
            "Length": len(read1.seq),
            "Name": "Unmapped"}
        list_names.append("Unmapped")
        idx_cigar += 1
    else:
        for cigar_tuple in read1.cigartuples:
            ad_name = CIGAR_DEFS[cigar_tuple[0]]
            if ad_name in list_names:
                ad_name = ad_name + ".2"
            else:
                list_names.append(ad_name)
            cigar_stat[idx_cigar] = {
                "Name": ad_name,
                "Length": cigar_tuple[1]}
            idx_cigar += 1
    if read2.cigarstring is None:
        cigar_stat[idx_cigar] = {
            "Length": len(read2.seq),
            "Name": "Unmapped"}
        idx_cigar += 1
    else:
        for cigar_tuple in read2.cigartuples:
            ad_name = CIGAR_DEFS[cigar_tuple[0]]
            if ad_name in list_names:
                ad_name = ad_name + ".2"
            else:
                list_names.append(ad_name)
            cigar_stat[idx_cigar] = {
                "Name": ad_name,
                "Length": cigar_tuple[1]}
            idx_cigar += 1
    return cigar_stat


def check_secondary_and_cigar(read1, read2):
    states = []
    for each_read in [read1, read2]:
        softclip = False
        for each_tuple in each_read.cigar:
            if each_tuple[0] == 4:
                softclip = True
        if each_read.is_secondary and softclip:
            states.append(False)
        elif softclip and not each_read.is_secondary:
            states.append(True)
        elif not each_read.is_secondary:
            states.append(False)
        else:
            states.append(False)
    if True in states:
        return True
    else:
        return False


def get_read_dict_paired(bam_path, other_bam_path="NA"):
    print("Working on {}".format(bam_path))
    sam_link = pysam.AlignmentFile(bam_path, "rb")
    # order_chroms = []
    dict_reads = {}
    read1 = None
    read2 = None
    read1name = ""
    read2name = " "
    reads_no_sequence = []
    for each_read in sam_link:
        if each_read.is_read2:
            read2 = each_read
            read2name = read2.qname
        else:
            read1 = each_read
            read1name = read1.qname
            read2 = None
        same_read = read1name == read2name
        if read1 is not None and read2 is not None and same_read:
            keep_read = check_secondary_and_cigar(read1, read2)
            if keep_read:
                idx_read = 1
                for read in [read1, read2]:
                    read_name = "{}.{}".format(read1.qname, idx_read)
                    dict_reads = parse_read(
                        read, dict_reads, read_name)
                    idx_read = idx_read + 1
    for each_read, each_dict in dict_reads.items():
        if each_dict.get("Sequence", None) is None:
            reads_no_sequence.append(each_read)
    sam_link = pysam.AlignmentFile(bam_path, "rb")
    for each_read in sam_link:
        if each_read.qname in reads_no_sequence:
            temp_dict = parse_read(
                each_read, {}, each_read.qname)
            cur_seq = temp_dict[each_read.qname].get("Sequence", None)
            if cur_seq is not None:
                dict_reads[each_read.qname]["Sequence"] = cur_seq
    return dict_reads


def score_motor(dict_virus, cigar_virus, len_seq):
    score_virus = np.zeros(len_seq)
    st = 0
    for i in range(len(cigar_virus.keys())):
        end = cigar_virus[i]["Length"]
        if cigar_virus[i]["Name"] == "Match":
            if dict_virus['Negative strand']:
                if st != 0:
                    score_virus[-end:-st] = 1
                else:
                    score_virus[-end:] = 1
            else:
                score_virus[st:end] = 1
        st = st + end
    return score_virus


def get_read_score(dict_virus, dict_host):
    cigar_virus = dict_virus["CIGAR"]
    len_seq = 0
    for i in range(len(cigar_virus.keys())):
        len_seq = len_seq + cigar_virus[i]["Length"]
    score_virus = score_motor(dict_virus, cigar_virus, len_seq)
    cigar_host = dict_host["CIGAR"]
    score_host = score_motor(dict_host, cigar_host, len_seq)
    xor_ar = np.logical_xor(score_virus, score_host)
    and_ar = np.logical_and(score_virus, score_host)
    out_score = (sum(xor_ar) - sum(and_ar)) / xor_ar.shape
    return out_score


def get_pos_strand(dict_virus):
    chrom = dict_virus["Chrom"]
    if dict_virus["Negative strand"]:
        strand = "Negative"
        start = dict_virus["Start"]
    else:
        strand = "Positive"
        start = dict_virus["End"]
    return chrom, start, strand


def get_ar_paired(dict_read_1, dict_read_2):
    '''
    Generates an array with paired-end information
    '''
    if dict_read_1.get("Sequence", None) is None:
        dict_read_1["Sequence"] = ""
    if dict_read_2.get("Sequence", None) is None:
        dict_read_2["Sequence"] = ""
    min_start = min(dict_read_1.get("Start", np.inf),
                    dict_read_2.get("Start", np.inf))
    max_end = max(dict_read_1.get("End", -np.inf),
                  dict_read_2.get("End", -np.inf))
    seq_length_1 = len(dict_read_1.get("Sequence", ""))
    seq_length_2 = len(dict_read_2.get("Sequence", ""))
    np_ar_temp = np.zeros(max_end + seq_length_1 + seq_length_2)
    np_ar_seq = np.empty(np_ar_temp.shape[0], dtype="|S16")
    last_reg = 0
    inverse = False
    strand = "Positive"
    chrom = np.nan
    for each_dict in [dict_read_1, dict_read_2]:
        new_chrom = each_dict.get("Chrom", np.nan)
        if "Chrom" in each_dict.keys():
            chrom = new_chrom
        start = each_dict.get("Start", np.nan)
        cigar_dict = each_dict.get("CIGAR", {})
        cur_st = start
        cur_seq = each_dict.get("Sequence", "")
        other_seq = each_dict.get("OtherSeq", "")
        if each_dict.get("Negative strand", False):
            strand = "Negative"
        if other_seq is not None and len(other_seq) > 0:
            _, rev_comp = make_rev_comp(cur_seq)
            if rev_comp == other_seq:
                inverse = True
        idx_seq = 0
        for i in range(len(cigar_dict)):
            len_reg = cigar_dict[i]["Length"]
            name_reg = cigar_dict[i]["Name"]
            ad_val = 0
            if "Soft" in name_reg:
                ad_val = -1
            if "Match" in name_reg:
                ad_val = 1
            if ad_val in [1, -1]:
                np_ar_temp[cur_st:cur_st + len_reg] = ad_val
                ad_seq = list(
                    cur_seq[idx_seq:idx_seq + len_reg])
                np_ar_seq[cur_st:cur_st + len(ad_seq)] = ad_seq
                idx_seq = idx_seq + len_reg
            last_reg = max(last_reg, cur_st + len_reg)
            cur_st = cur_st + len_reg
    np_ar_temp = np_ar_temp[min_start:last_reg]
    np_ar_seq = np_ar_seq[min_start:last_reg]
    if inverse:
        np_ar_temp = np_ar_temp[::-1]
        np_ar_seq = np_ar_seq[::-1]
    return np_ar_temp, min_start, np_ar_seq, strand, chrom


def match_sequences(dict_virus_1, dict_host_1):
    if dict_virus_1.get("Sequence", None) is None:
        strand1 = dict_virus_1.get("Negative strand", -1)
        if dict_host_1.get("Sequence", None) is not None:
            ad_seq = dict_host_1["Sequence"]
            strand2 = dict_host_1.get("Negative strand", -1)
            if np.logical_and(strand1, strand2):
                dict_virus_1["Sequence"] = ad_seq
            elif strand2:
                _, ad_seq_rc = make_rev_comp(ad_seq)
                dict_virus_1["Sequence"] = ad_seq_rc
    elif dict_host_1.get("Sequence", None) is None:
        strand1 = dict_host_1.get("Negative strand", -1)
        if dict_virus_1.get("Sequence", None) is not None:
            ad_seq = dict_virus_1["Sequence"]
            strand2 = dict_virus_1.get("Negative strand", -1)
            if np.logical_and(strand1, strand2):
                dict_host_1["Sequence"] = ad_seq
            elif strand2:
                _, ad_seq_rc = make_rev_comp(ad_seq)
                dict_host_1["Sequence"] = ad_seq_rc
    return dict_virus_1, dict_host_1


def get_paired_score(dict_virus_1, dict_virus_2, dict_host_1, dict_host_2):
    dict_virus_1, dict_host_1 = match_sequences(dict_virus_1, dict_host_1)
    dict_virus_2, dict_host_2 = match_sequences(dict_virus_2, dict_host_2)
    dict_virus_1["OtherSeq"] = dict_host_1.get("Sequence", "")
    dict_virus_2["OtherSeq"] = dict_host_2.get("Sequence", "")
    paired_ar_vir, start_vir, seq_vir, strand_vir, chr_v = get_ar_paired(
        dict_virus_1, dict_virus_2)
    paired_ar_host, start_host, seq_host, strand_host, chr_h = get_ar_paired(
        dict_host_1, dict_host_2)
    # Will only the information for the region un/mapped by both (0 to min_len)
    # If a paired read is not mapped to one, will therefore be excluded
    min_len = min(len(paired_ar_vir), len(paired_ar_host))
    out_ar = np.logical_not(paired_ar_vir[:min_len] + paired_ar_host[:min_len])
    out_score = np.mean(out_ar)
    read_name = dict_host_1.get("ReadName", "")
    if read_name == "":
        read_name = dict_host_2.get("ReadName", "")
    dict_out = {"Score": out_score,
                "Start.Host": start_host,
                "Start.Virus": start_vir,
                "Strand.Virus": strand_vir,
                "Strand.Host": strand_host,
                "Chrom.Host": chr_h,
                "Chrom.Virus": chr_v,
                "ReadName": read_name}
    return dict_out


def annotate_regions_paired(dict_reads_virus, dict_reads_host):
    dict_integration = {}
    for each_read in dict_reads_virus.keys():
        if each_read.split(".")[-1] == "1":
            read_1 = each_read
            read_2 = ".".join(each_read.split(".")[:-1]) + ".2"
        else:
            read_2 = each_read
            read_1 = ".".join(each_read.split(".")[:-1]) + ".1"
        found_read_1 = read_1 in dict_reads_host.keys()
        found_read_2 = read_2 in dict_reads_host.keys()
        if found_read_1 or found_read_2:
            # print("Investigating {}".format(read_1))
            dict_virus_1 = dict_reads_virus.get(read_1, {})
            dict_virus_2 = dict_reads_virus.get(read_2, {})
            not_second_virus = np.logical_or(
                dict_virus_1.get("IsPrimary", False),
                dict_virus_2.get("IsPrimary", False))
            has_vir_softclip = np.logical_or(
                len(dict_virus_1.get("CIGAR", {})) > 1,
                len(dict_virus_2.get("CIGAR", {})) > 1)
            dict_host_1 = dict_reads_host.get(read_1, {})
            dict_host_2 = dict_reads_host.get(read_2, {})
            not_second_host = np.logical_or(
                dict_host_1.get("IsPrimary", False),
                dict_host_2.get("IsPrimary", False))
            has_host_softclip = np.logical_or(
                len(dict_host_1.get("CIGAR", {})) > 1,
                len(dict_host_2.get("CIGAR", {})) > 1)
            has_2_softclip = np.logical_and(
                len(dict_host_2.get("CIGAR", {})) > 1,
                len(dict_virus_2.get("CIGAR", {})) > 1)
            has_1_softclip = np.logical_and(
                len(dict_host_1.get("CIGAR", {})) > 1,
                len(dict_virus_1.get("CIGAR", {})) > 1)
            has_match_softclip = has_1_softclip or has_2_softclip
            has_both_softclip = has_vir_softclip and has_host_softclip
            has_primary = np.logical_and(
                not_second_virus, not_second_host)
            if has_match_softclip and has_both_softclip:
                if has_primary:
                    dict_integration[each_read] = get_paired_score(
                        dict_virus_1, dict_virus_2,
                        dict_host_1, dict_host_2)
    return dict_integration


def annotate_regions(dict_reads_virus, dict_reads_host):
    dict_integration = {}
    for each_read in dict_reads_virus.keys():
        if each_read in dict_reads_host.keys():
            dict_virus = dict_reads_virus[each_read]
            dict_host = dict_reads_host[each_read]
            score_read = get_read_score(dict_virus, dict_host)
            chr_vr, pos_vr, strand_vr = get_pos_strand(
                dict_virus)
            chr_h, pos_h, strand_h = get_pos_strand(
                dict_host)
            dict_integration[each_read] = {
                "Chrom.Host": chr_h,
                "Start.Host": pos_h,
                "Strand.Host": strand_h,
                "Chrom.Virus": chr_vr,
                "Start.Virus": pos_vr,
                "Strand.Virus": strand_vr,
                "Score": score_read,
                "ReadName": each_read}
    return dict_integration


def get_chrom_dict_integ(integration_dict, window):
    chroms = {}
    idx_read = 0
    for each_read, each_dict in integration_dict.items():
        chrom = each_dict["Chrom.Host"]
        pos = each_dict["Start.Host"]
        if chrom not in chroms.keys():
            range_pos = range(pos - window, pos + window)
            chroms[chrom] = {"Poses": [pos],
                             "Range": range_pos}
        else:
            if pos not in chroms[chrom]["Range"]:
                cur_list = chroms[chrom]["Poses"]
                cur_list.append(pos)
                chroms[chrom]["Poses"] = cur_list
        idx_read += 1
    return chroms


def summarize_integration(integration_dict, window):
    used_reads = []
    summary_dict = {}
    chroms = get_chrom_dict_integ(integration_dict, window)
    for chrom in chroms.keys():
        print("Working on {} with {} integrations".format(
            chrom, len(chroms[chrom]["Poses"])))
        for each_pos in chroms[chrom]["Poses"]:
            range_pos = range(each_pos - window, each_pos + window)
            for each_read, each_dict in integration_dict.items():
                cur_pos = each_dict["Start.Host"]
                cur_chrom = each_dict["Chrom.Host"]
                if cur_chrom == chrom and cur_pos in range_pos:
                    if each_read not in used_reads:
                        used_reads.append(each_read)
                        ad_name = "{}.{}".format(chrom, each_pos)
                        cur_dict = summary_dict.get(ad_name, {})
                        if "Chrom.Host" not in cur_dict.keys():
                            ad_names = ['Chrom.Host', 'Start.Host',
                                        'Strand.Host', 'Chrom.Virus',
                                        'Start.Virus', 'Strand.Virus',
                                        'ReadName']
                            for each in ad_names:
                                cur_dict[each] = str(each_dict[each])
                            cur_dict["Supporting reads"] = 1
                            cur_dict["Score"] = each_dict["Score"]
                        else:
                            ad_names = ['Start.Host', 'Strand.Host',
                                        'Start.Virus', 'Strand.Virus',
                                        'ReadName']
                            for each in ad_names:
                                cur_dict[each] = "{}, {}".format(
                                    cur_dict[each], each_dict[each])
                            cur_dict["Supporting reads"] += 1
                            cur_dict["Score"] = cur_dict["Score"] +\
                                each_dict["Score"]
                        summary_dict[ad_name] = cur_dict
    return summary_dict


def find_integration(host_bam, viral_bam, paired_end):
    WINDOW = 25
    print("Started working on {} and {}".format(host_bam, viral_bam))
    if paired_end:
        dict_reads_virus = get_read_dict_paired(viral_bam, host_bam)
        dict_reads_host = get_read_dict_paired(host_bam, viral_bam)
        integration_dict = annotate_regions_paired(
            dict_reads_virus, dict_reads_host)
    else:
        dict_reads_virus = get_read_dict(viral_bam)
        dict_reads_host = get_read_dict(host_bam)
        integration_dict = annotate_regions(
            dict_reads_virus, dict_reads_host)
    print(integration_dict)
    out_dict = summarize_integration(
        integration_dict, WINDOW)
    out_list = []
    for each_integration, each_dict in out_dict.items():
        if int(each_dict["Supporting reads"]) > 1:
            read_names = each_dict.get("ReadName", "NoReadFound")
            read_names = ", ".join(np.unique(read_names.split(", ")))
            supp_reads = each_dict["Supporting reads"]
            score_reads = each_dict["Score"]
            try:
                supp_reads = supp_reads[0]
            except Exception:
                None
            try:
                score_reads = score_reads[0]
            except Exception:
                None
            ad_list = [
                host_bam, viral_bam, each_dict["Chrom.Host"],
                each_dict["Start.Host"], each_dict["Strand.Host"],
                each_dict["Chrom.Virus"], each_dict["Start.Virus"],
                each_dict["Strand.Virus"],
                str(supp_reads),
                str(score_reads),
                str(score_reads / supp_reads),
                read_names]
            out_list.append(ad_list)
    return out_list


class polyidusEngine:
    def __init__(self, hostindex, viralindex,
                 fastq, outdir, aligner, virname):
        self.hostindex = hostindex
        self.viralindex = viralindex
        self.fastq = fastq
        self.paired = False
        if len(self.fastq) == 2:
            self.paired = True
        self.outdir = outdir
        self.aligner = aligner
        self.virname = virname
        self.make_folders()
        self.viralbam_temp = os.path.join(
            self.outdir_viral, "virusAligned-temp.bam")
        self.viralbam_final = os.path.join(
            self.outdir_viral, "virusAligned-final.bam")
        self.hostbam = os.path.join(self.outdir_host, "hostAligned.bam")
        self.hostbam_sorted = os.path.join(
            self.outdir_host, "hostAligned_sorted.bam")
        self.command_lists = []

    def update_log(self):
        logpath = os.path.join(
            self.outdir,
            "polyidusCommands.sh")
        with open(logpath, "w") as loglink:
            for each in self.command_lists:
                outstr = " ".join([str(val) for val in each]) + "\n"
                loglink.write(outstr)

    def make_folders(self):
        fold_dict = {"outdir_viral": "viral",
                     "outdir_host": "host",
                     "outdir_results": "results"}
        for varname, foldername in fold_dict.items():
            outdir = os.path.join(self.outdir, foldername)
            os.makedirs(outdir, exist_ok=True)
            setattr(self, varname, outdir)

    def align_files(self):
        print(
            "Aligning to viral genome at: {}".format(
                str(datetime.now())))
        if self.aligner == "bwa":
            self.align_virus_bwa()
        else:
            self.align_virus()
        print(
            "Extracting mapped fastqs from the viral genome at: {}".format(
                str(datetime.now())))
        self.extract_virus()
        print(
            "Aligning viral-mapped FASTQs to host at: {}".format(
                str(datetime.now())))
        if self.aligner == "bwa":
            self.align_host_bwa()
        else:
            self.align_host()
        print(
            "All BAM files generated successfully at: {}".format(
                str(datetime.now())))

    def align_virus(self):
        if len(self.fastq) == 1:
            job1 = [
                "bowtie2", "-p", "1", "--local", "-x", self.viralindex,
                "-U", self.fastq[0]]
            job2 = [
                "samtools",
                "view", "-bS", "-F", "4", "-o",
                self.viralbam_final,
                "-"]
        elif len(self.fastq) == 2:
            job1 = [
                "bowtie2", "-p", "1", "--local",
                "-x", self.viralindex, "-1",
                self.fastq[0], "-2", self.fastq[1]]
            job2 = [
                "samtools", "sort", "-o",
                self.viralbam_temp, "-"]
        self.command_lists.append(job1 + ["|"] + job2)
        self.update_log()
        p1 = subprocess.Popen(job1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(job2, stdin=p1.stdout, stdout=subprocess.PIPE)
        p2.communicate()
        if len(self.fastq) == 2:
            viralbam_temp1 = os.path.join(
                self.outdir_viral, "pair1mapped.bam")
            viralbam_temp2 = os.path.join(
                self.outdir_viral, "pair2mapped.bam")
            viralbam_bothmapped = os.path.join(
                self.outdir_viral, "bothPairsMapped.bam")
            # Part 1. Pair 1 mapped pair 2 unmapped
            job3 = [
                "samtools", "view", "-bS", "-f",
                "4", "-F", "264", self.viralbam_temp,
                "-o", viralbam_temp1]
            # Part 2. Pair 2 mapped pair 1 unmapped
            job4 = [
                "samtools", "view", "-bS", "-f", "8", "-F", "260",
                self.viralbam_temp, "-o", viralbam_temp2]
            # Pair 3. both pairs mapped
            job5 = [
                "samtools", "view", "-bS", "-f", "1", "-F", "12",
                self.viralbam_temp, "-o", viralbam_bothmapped]
            for eachjob in [job3, job4, job5]:
                self.command_lists.append(eachjob)
                self.update_log()
                subprocess.run(eachjob, check=True)
            job6 = ["samtools", "merge", "-", viralbam_temp1,
                    viralbam_temp2, viralbam_bothmapped]
            job7 = ["samtools", "sort", "-n", "-", "-o", self.viralbam_final]
            self.command_lists.append(job6 + ["|"] + job7)
            self.update_log()
            p3 = subprocess.Popen(job6, stdout=subprocess.PIPE)
            p4 = subprocess.Popen(
                job7, stdin=p3.stdout, stdout=subprocess.PIPE)
            p4.communicate()

    def extract_virus(self):
        self.fastq_path_1 = os.path.join(
            self.outdir_viral, "ViralAligned_1.fastq")
        job1 = ["bedtools", "bamtofastq", "-i",
                self.viralbam_final, "-fq", self.fastq_path_1]
        if len(self.fastq) == 2:
            self.fastq_path_2 = os.path.join(
                self.outdir_viral, "ViralAligned_2.fastq")
            job1.extend(["-fq2", self.fastq_path_2])
        self.command_lists.append(job1)
        self.update_log()
        subprocess.run(job1, check=True)

    def align_host(self):
        if len(self.fastq) == 1:
            job1 = [
                "bowtie2", "-p", "1", "--local", "-x", self.hostindex,
                "U", self.fastq_path_1]
            job2 = [
                "samtools",
                "view", "-bS", "-F", "4", "-o",
                self.hostbam,
                "-"]
        elif len(self.fastq) == 2:
            job1 = [
                "bowtie2", "-p", "1", "--local",
                "-x", self.hostindex, "-1",
                self.fastq_path_1, "-2", self.fastq_path_2]
            # Version 1.1.0: sort by name! default was position :(
            job2 = [
                "samtools", "sort", "-n", "-o",
                self.hostbam, "-"]
        self.command_lists.append(job1 + ["|"] + job2)
        self.update_log()
        p1 = subprocess.Popen(job1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(job2, stdin=p1.stdout, stdout=subprocess.PIPE)
        p2.communicate()
        # Sort host output
        job3 = ["samtools", "sort", self.hostbam, "-o", self.hostbam_sorted]
        job4 = ["samtools", "index", self.hostbam_sorted]
        for eachjob in [job3, job4]:
            self.command_lists.append(eachjob)
            self.update_log()
            subprocess.run(eachjob)

    def align_virus_bwa(self):
        if len(self.fastq) == 1:
            job1 = [
                "bwa", "mem", "-t", "1", "-T", "1", "-a",
                "-C", "-Y", self.viralindex,
                self.fastq[0]]
            job2 = [
                "samtools",
                "view", "-bS", "-F", "4", "-o",
                self.viralbam_final,
                "-"]
        elif len(self.fastq) == 2:
            job1 = [
                "bwa", "mem", "-t", "1", "-T", "1", "-a", "-C", "-Y",
                self.viralindex,
                self.fastq[0], self.fastq[1]]
            job2 = [
                "samtools", "sort", "-o",
                self.viralbam_temp, "-"]
        self.command_lists.append(job1 + ["|"] + job2)
        self.update_log()
        p1 = subprocess.Popen(job1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(job2, stdin=p1.stdout, stdout=subprocess.PIPE)
        p2.communicate()
        if len(self.fastq) == 2:
            viralbam_temp1 = os.path.join(
                self.outdir_viral, "pair1mapped.bam")
            viralbam_temp2 = os.path.join(
                self.outdir_viral, "pair2mapped.bam")
            viralbam_bothmapped = os.path.join(
                self.outdir_viral, "bothPairsMapped.bam")
            job3 = [
                "samtools", "view", "-bS", "-f",
                "4", "-F", "264", self.viralbam_temp,
                "-o", viralbam_temp1]
            job4 = [
                "samtools" "view", "-bS", "-f", "8", "-F", "260",
                self.viralbam_temp, "-o", viralbam_temp2]
            job5 = [
                "samtools", "view", "-bS", "-f", "1", "-F", "12",
                self.viralbam_temp, "-o", viralbam_bothmapped]
            for eachjob in [job3, job4, job5]:
                self.command_lists.append(eachjob)
                self.update_log()
                subprocess.run(eachjob, check=True)
            job6 = ["samtools", "merge", "-", viralbam_temp1,
                    viralbam_temp2, viralbam_bothmapped]
            job7 = ["samtools", "sort", "-n", "-", "-o", self.viralbam_final]
            self.command_lists.append(job6 + ["|"] + job7)
            self.update_log()
            p3 = subprocess.Popen(job6, stdout=subprocess.PIPE)
            p4 = subprocess.Popen(
                job7, stdin=p3.stdout, stdout=subprocess.PIPE)
            p4.communicate()

    def align_host_bwa(self):
        if len(self.fastq) == 1:
            job1 = [
                "bwa", "mem", "-t", "4", "-T",
                "1", "-a", "-C", "-Y", self.hostindex,
                self.fastq_path_1]
            job2 = [
                "samtools",
                "view", "-bS", "-F", "4", "-o",
                self.hostbam,
                "-"]
        elif len(self.fastq) == 2:
            job1 = [
                "bwa", "mem", "-t", "4", "-T",
                "1", "-a", "-C", "-Y",
                self.hostindex,
                self.fastq_path_1, self.fastq_path_2]
            job2 = [
                "samtools", "sort", "-o",
                self.hostbam, "-"]
        self.command_lists.append(job1 + ["|"] + job2)
        self.update_log()
        p1 = subprocess.Popen(job1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(job2, stdin=p1.stdout, stdout=subprocess.PIPE)
        p2.communicate()
        # Sort host output
        job3 = ["samtools", "sort", self.hostbam, "-o", self.hostbam_sorted]
        job4 = ["samtools", "index", self.hostbam_sorted]
        for eachjob in [job3, job4]:
            self.command_lists.append(eachjob)
            self.update_log()
            subprocess.run(eachjob)

    def find_approximate_integrations(self):
        self.integpath_approximate = os.path.join(
            self.outdir_results,
            "{}IntegrationInfo.tsv".format(self.virname))
        integlist = [["HostFile", "ViralFile", "ChromHost",
                      "PositionHost", "StrandHost",
                      "ChromViral", "PositionViral", "StrandViral",
                      "NumberReads", "Score",
                      "Normalized.Score", "ReadNames"]]
        ad_list = find_integration(
            self.hostbam, self.viralbam_final, self.paired)
        integlist.extend(ad_list)
        with open(self.integpath_approximate, "w") as out_link:
            for each_list in integlist:
                out_str = "\t".join(each_list) + "\n"
                print(out_str)
                out_link.write(out_str)

    def find_exact_integrations(self):
        # Defaults
        minreads = 2
        minscore = 0.6
        # Get sample name
        path = os.path.normpath(self.viralbam_final)
        pathparts = path.split(os.sep)
        samplename = pathparts[-3]
        # Make output path
        integpath = os.path.join(
            self.outdir_results,
            "exact{}Integrations.tsv".format(self.virname))
        integpathbed = os.path.join(
            self.outdir_results,
            "exact{}Integrations.bed".format(self.virname))
        outlink = open(integpath, "w")
        bedlink = open(integpathbed, "w")
        headerbed = ["chrom", "chromStart", "chromEnd",
                     "name", "score", "strand",
                     "virusStart", "virusEnd",
                     "virusStrand", "numReads", "bioSample"]
        bedlink.write("\t".join(headerbed) + "\n")
        header = ["Chrom", "Start", "End", "IntegrationSite",
                  "ChromVirus", "VirusStart", "VirusEnd",
                  "ViralIntegrationSite",
                  "NumberOfSupportingFragments", "SampleName",
                  "FragmentName"]
        keys = ["ChromHost", "Starts", "Ends", "IntegrationHosts",
                "ChromVirus", "VirusStarts", "VirusEnds",
                "IntegrationVirus", "NumberReads"]
        outlink.write("\t".join(header) + "\n")
        chimerdf = get_chimer(self.integpath_approximate, minreads, minscore)
        for i in range(chimerdf.shape[0]):
            chromhost = str(chimerdf.iloc[i, 2])
            if "alt" not in chromhost and "rand" not in chromhost:
                readnames = chimerdf.iloc[i, -1].split(", ")
                start_host = int(
                    np.mean(
                        [int(each) for each in
                         chimerdf.iloc[i, 3].split(", ")]))
                dictposes = get_host_pos(
                    self.viralbam_final, self.hostbam_sorted, chromhost,
                    start_host, readnames)
                LEN_VIR = len(dictposes["ChromVirus"]) > 0
                if LEN_VIR and dictposes["ChromVirus"][0] != "NA":
                    for j in range(len(dictposes["Starts"])):
                        adlist = []
                        for key in keys:
                            adlist.append(str(dictposes[key][j]))
                        adlist.append(samplename)
                        fragnames = dictposes[
                            "ReadNameDict"][dictposes["IntegrationHosts"][j]]
                        adlist.append(", ".join(fragnames))
                        print("\t".join(adlist))
                        outlink.write("\t".join(adlist) + "\n")
                        strand_host = "+"
                        strand_virus = "+"
                        if dictposes["Starts"][j] > \
                                dictposes["IntegrationHosts"][j]:
                            strand_host = "-"
                        if dictposes["VirusStarts"][j] > \
                                dictposes["IntegrationVirus"][j]:
                            strand_virus = "-"
                        adlistbed = [dictposes["ChromHost"][j],
                                     dictposes["IntegrationHosts"][j] - 1,
                                     dictposes["IntegrationHosts"][j],
                                     dictposes["ChromVirus"][j], "0",
                                     strand_host,
                                     dictposes["IntegrationVirus"][j] - 1,
                                     dictposes["IntegrationVirus"][j],
                                     strand_virus,
                                     dictposes["NumberReads"][j], samplename]
                        adlistbed = [str(each) for each in adlistbed]
                        bedlink.write("\t".join(adlistbed) + "\n")
        outlink.close()
        bedlink.close()
