#!/usr/bin/python
#-*- coding: UTF-8 -*-

import os
import sys
from optparse import OptionParser
import pysam
import time
import multiprocessing


def parseCommand():
    usage = "usage: ./dedup.py -1 in.sorted.bam -o out.bam > dedup.log"
    version = "%prog 1.0"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-1", "--input1", dest="input1", help="the bam file")
    parser.add_option("-q", "--phredQ", dest="phred",
                      help="Phred quality score threshold", type="int", default="20")
    parser.add_option("-o", "--output", dest="output",
                      default=".", help="output file".decode("utf8"))
    return parser.parse_args()


def str_compare(forward_list1, forward_list2, reverse_list1, reverse_list2, phredQ=20):
    seq_1 = forward_list1[2] + reverse_list1[2]
    phred_1 = forward_list1[3].tolist()
    tmp_2 = reverse_list1[3].tolist()
    phred_1.extend(tmp_2)
    seq_2 = forward_list2[2] + reverse_list2[2]
    phred_2 = forward_list2[3].tolist()
    tmp_2 = reverse_list2[3].tolist()
    phred_2.extend(tmp_2)
    if len(seq_1) == len(seq_2):
        if seq_1 == seq_2:
            return True
        else:
            index = [0]
            count = 0
            for i in range(len(seq_1)):
                if phred_1[i] < phredQ or phred_2[i] < phredQ:
                    if i == 0:
                        if seq_1[i] != seq_2[i]:
                            return False
                    else:
                        index.append(i)
                        count += 1
                        if seq_1[(index[count - 1] + 1):index[count]] != seq_2[(index[count - 1] + 1):index[count]]:
                            return False
                else:
                    if seq_1[i] == seq_2[i]:
                        continue
                    else:
                        return False
            return True
    else:
        return False


def block_compare(forward_block_dict_list, reverse_block_dict, phredQ=20):
    dup_seqname_dict = {}
    save_number_dict = {}
    line_number_list = [save_number_dict, dup_seqname_dict]
    start = True
    for i in range(len(forward_block_dict_list)):
        if forward_block_dict_list[i][1] not in dup_seqname_dict:
            save_number_dict[forward_block_dict_list[i][0]] = True
            if forward_block_dict_list[i][1] not in reverse_block_dict:
                continue
        for j in range((i + 1), len(forward_block_dict_list)):
            if forward_block_dict_list[j][1] not in dup_seqname_dict:
                if forward_block_dict_list[j][1] in reverse_block_dict:
                    if str_compare(forward_block_dict_list[i], forward_block_dict_list[j], reverse_block_dict[forward_block_dict_list[i][1]], reverse_block_dict[forward_block_dict_list[j][1]], phredQ):
                        dup_seqname_dict[forward_block_dict_list[j][1]] = True
    return line_number_list


def dedup3(in_file, phred, out_file):
    f = pysam.AlignmentFile(in_file, "rb")
    total = 0
    buffer_block = []
    forward_block_dict = {}
    reverse_block_dict = {}
    other_dict = {}
    chrom_flag = False
    chrom_tmp = "chrM"
    line_list = []
    dedup_dict = {}
    save_dict = {}
    for frag in f:
        if frag.is_proper_pair:
            if frag.template_length > 0:
                if frag.reference_id not in forward_block_dict:
                    forward_block_dict[frag.reference_id] = {}
                    forward_block_dict[frag.reference_id][frag.reference_start] = [
                        [total, frag.query_name, frag.query_sequence, frag.query_qualities]]
                    if not chrom_flag:
                        chrom_tmp = frag.reference_id
                        chrom_flag = True
                        total += 1
                        continue
                    for i in forward_block_dict[chrom_tmp].keys():
                        if len(forward_block_dict[chrom_tmp][i]) >= 2:
                            line_list = block_compare(forward_block_dict[chrom_tmp][
                                                      i], reverse_block_dict[chrom_tmp], phred)
                            dedup_dict.update(line_list[1])
                            save_dict.update(line_list[0])
                        else:
                            save_dict[forward_block_dict[
                                chrom_tmp][i][0][0]] = True
                    forward_block_dict[chrom_tmp].clear()
                    for i in reverse_block_dict[chrom_tmp].keys():
                        if i not in dedup_dict:
                            save_dict[reverse_block_dict[
                                chrom_tmp][i][0]] = True
                    reverse_block_dict[chrom_tmp].clear()
                    chrom_tmp = frag.reference_id
                elif frag.reference_start not in forward_block_dict[frag.reference_id]:
                    forward_block_dict[frag.reference_id][
                        frag.reference_start] = {}
                    forward_block_dict[frag.reference_id][frag.reference_start] = [
                        [total, frag.query_name, frag.query_sequence, frag.query_qualities]]
                else:
                    forward_block_dict[frag.reference_id][frag.reference_start].append(
                        [total, frag.query_name, frag.query_sequence, frag.query_qualities])
            elif frag.template_length < 0:
                if frag.reference_id not in reverse_block_dict:
                    reverse_block_dict[frag.reference_id] = {}
                    reverse_block_dict[frag.reference_id][frag.query_name] = [
                        total, frag.query_name, frag.query_sequence, frag.query_qualities]
                else:
                    reverse_block_dict[frag.reference_id][frag.query_name] = [
                        total, frag.query_name, frag.query_sequence, frag.query_qualities]
            else:
                other_dict[total] = True
        else:
            other_dict[total] = True
        total += 1

    if chrom_tmp in forward_block_dict:
        for i in forward_block_dict[chrom_tmp].keys():
            if len(forward_block_dict[chrom_tmp][i]) >= 2:
                line_list = block_compare(forward_block_dict[chrom_tmp][
                                          i], reverse_block_dict[chrom_tmp], phred)
                dedup_dict.update(line_list[1])
                save_dict.update(line_list[0])
            else:
                save_dict[forward_block_dict[chrom_tmp][i][0][0]] = True
        forward_block_dict[chrom_tmp].clear()
        for i in reverse_block_dict[chrom_tmp].keys():
            if i not in dedup_dict:
                save_dict[reverse_block_dict[chrom_tmp][i][0]] = True
        reverse_block_dict[chrom_tmp].clear()

    save_dict.update(other_dict)
    f.close()
    f = pysam.AlignmentFile(in_file, "rb")
    f_w = pysam.AlignmentFile(out_file, "wb", template=f)
    print "total=", total,  "---after dedup---", len(save_dict.keys())
    total = 0
    for frag in f:
        if total in save_dict:
            f_w.write(frag)
        total += 1
    f_w.close()
    print "total=", total,  "---after dedup---", len(save_dict.keys())


if __name__ == "__main__":
    (options, args) = parseCommand()
    if options.input1 == None:
        print "see -h for help"
        sys.exit(-1)
    if options.output == None:
        print "see -h for help"
        sys.exit(-1)
    time1 = time.time()
    dedup3(options.input1, options.phred, options.output)
    time2 = time.time()
    print 'Time used: ' + str(time2 - time1)