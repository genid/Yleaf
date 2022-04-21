#!/usr/bin/env python

"""
Yleaf detection of Y-Haplogroups in Human DNA v3.0

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Diego Montiel Gonzalez
Slightly modified by: Bram van Wersch
"""

import os
import sys
import re
import time
import subprocess
import pandas as pd
import numpy as np
from argparse import ArgumentParser
import gc
from tree import Tree

pd.options.mode.chained_assignment = None  # default='warn'

VERSION = 3.0


def get_arguments():
    parser = ArgumentParser()

    parser.add_argument("-fastq", "--fastq",
                        dest="Fastq", required=False,
                        help="Use raw FastQ files", metavar="PATH")

    parser.add_argument("-bam", "--bam",
                        dest="Bamfile", required=False,
                        help="input BAM file", metavar="PATH")

    parser.add_argument("-cram", "--cram",
                        dest="Cramfile", required=False,
                        help="input CRAM file", metavar="PATH")

    parser.add_argument("-f", "--fasta-ref", dest="reference",
                        help="fasta reference genome sequence ", metavar="PATH", required=False)

    parser.add_argument("-pos", "--position", dest="position",
                        help="Positions file [hg19.txt or hg38.txt]", metavar="PATH", required=True)

    parser.add_argument("-out", "--output",
                        dest="Outputfile", required=True,
                        help="Folder name containing outputs", metavar="STRING")

    parser.add_argument("-r", "--Reads_thresh",
                        help="The minimum number of reads for each base",
                        type=int, required=False,
                        default=50)

    parser.add_argument("-q", "--Quality_thresh",
                        help="Minimum quality for each read, integer between 10 and 39, inclusive \n [10-40]",
                        type=int, required=True)

    parser.add_argument("-b", "--Base_majority",
                        help="The minimum percentage of a base result for acceptance \n [50-99]",
                        type=int, required=True)

    parser.add_argument("-t", "--Threads", dest="threads",
                        help="Set number of additional threads to use during alignment BWA-MEM",
                        type=int,
                        default=2)

    parser.add_argument("-old", "--use_old", dest="use_old",
                        help="Add this value if you want to use the old prediction method of Yleaf (version 2.3). This"
                             " version only uses the ISOGG tree and slightly different prediction criteria.",
                        action="store_true")

    args = parser.parse_args()
    return args


def get_frequency_table(mpileup):
    bases = ["A", "T", "G", "C", "+", "-"]
    frequency_table = {}

    for i in mpileup.values:
        fastadict = {"A": 0, "T": 0, "G": 0, "C": 0}
        sequence = i[9]  # actual sequence
        sequence = sequence.upper()
        sequence = trimm_caret(sequence)
        sequence = sequence.replace("$", "")
        indel_pos = find_all_indels(sequence)
        # Count number of indels
        indels = count_indels(sequence, indel_pos)
        fastadict.update(indels)
        fastadict["-"] += sequence.count("*")
        # Trimm Indels
        trimm_sequence = trimm_indels(sequence, indel_pos)
        for seq in trimm_sequence:
            if seq in fastadict:
                fastadict[seq] += 1
        frequency_table.update({i[3]: list(fastadict.values())})
    df_frequency_table = pd.DataFrame.from_dict(frequency_table, orient='index')
    df_frequency_table.columns = bases
    return df_frequency_table


def find_all_indels(s):
    list_pos = []
    for i in find_all(s, "-"):
        list_pos.append(i)
    for i in find_all(s, "+"):
        list_pos.append(i)
    return sorted(list_pos)


def count_indels(s, pos):
    dict_indel = {"+": 0, "-": 0}
    if len(pos) == 0:
        return dict_indel
    if len(pos) > 0:
        for i in range(0, len(pos)):
            try:  # in case it is not a number but a base pair e.g. A
                dict_indel[s[pos[i]]] += int(s[pos[i] + 1])
            except ValueError:
                dict_indel[s[pos[i]]] += 1
                continue
    return dict_indel


def trimm_indels(s, pos):
    # Receives a sequence and trimms indels
    if len(pos) == 0:
        return s
    start = pos[0]
    count = (start + 1)
    try:  # in case it is not a number but a base pair e.g. A
        end = count + int(s[count]) + 1
    except ValueError:
        end = start + 1
    u_sequence = s[:start]
    if len(pos) > 1:
        for i in range(1, len(pos)):
            start = end
            u_sequence += s[start:pos[i]]
            start = pos[i]
            count = (start + 1)
            try:  # in case it is not a number but a base pair e.g. A
                end = count + int(s[count]) + 1
            except ValueError:
                end = start + 1
            if pos[-1] == pos[i]:
                u_sequence += s[end:]
    else:
        u_sequence += s[end:]
    return u_sequence


def trimm_caret(s):
    list_pos = []
    for i in find_all(s, "^"):
        list_pos.append(i)
    if len(list_pos) == 0:
        return s
    i = 0
    start = 0
    sequence = ""
    while i < len(s):
        if s[i] == "^":
            end = i
            sequence += (s[start:end])
            start = i + 1
        elif i >= list_pos[-1] + 1:
            sequence += (s[list_pos[-1] + 1:])
            break
        i += 1
    return sequence


def find_all(c, s):
    return [x for x in range(c.find(s), len(c)) if c[x] == s]


def execute_mpileup(header, bam_file, pileupfile, quality_thresh, folder, reference):
    if reference:
        cmd = "samtools mpileup -f {} -AQ{} -r {} {} > {}".format(reference,
                                                                  quality_thresh, header, bam_file, pileupfile)

    else:
        cmd = "samtools mpileup -AQ{} -r {} {} > {}".format(quality_thresh,
                                                            header, bam_file, pileupfile)
    subprocess.call(cmd, shell=True)


def chromosome_table(path_file, path_folder, file_name):
    output = path_folder + '/' + file_name + '.chr'
    tmp_output = path_folder + "/tmp.txt"

    f = open(tmp_output, "w")
    subprocess.call(["samtools", "idxstats", path_file], stdout=f)
    df_chromosome = pd.read_table(tmp_output, header=None)

    total_reads = sum(df_chromosome[2])
    df_chromosome["perc"] = (df_chromosome[2] / total_reads) * 100
    df_chromosome = df_chromosome.round(decimals=2)
    df_chromosome['perc'] = df_chromosome['perc'].astype(str) + '%'
    df_chromosome = df_chromosome.drop(columns=[1, 3])
    df_chromosome.columns = ['chr', 'reads', 'perc']
    df_chromosome.to_csv(output, index=None, sep="\t")

    cmd = "rm " + tmp_output
    subprocess.call(cmd, shell=True)

    if 'Y' in df_chromosome["chr"].values:
        return "Y", total_reads
    elif 'chrY' in df_chromosome["chr"].values:
        return "chrY", total_reads


def check_if_folder(path, ext):
    list_files = []
    if os.path.isdir(path):
        dirpath = os.walk(path)
        for dirpath, dirnames, filenames in dirpath:
            for filename in [f for f in filenames if f.endswith(ext)]:
                files = os.path.join(dirpath, filename)
                list_files.append(files)
        return list_files
    else:
        return [path]


def get_folder_name(path_file):
    folder = path_file.split('/')[-1]
    folder_name = os.path.splitext(folder)[0]
    return folder_name


def create_tmp_dirs(folder):
    if os.path.isdir(folder):
        while True:
            print("WARNING! File " + folder + " already exists, \nWould you like to remove it?")
            choice = input("y/n: ")
            if str(choice).upper() == "Y":
                cmd = 'rm -r ' + folder
                subprocess.call(cmd, shell=True)
                cmd = 'mkdir ' + folder
                subprocess.call(cmd, shell=True)
                return True
            elif str(choice).upper() == "N":
                return False
            else:
                print("Please type y/Y or n/N")
    else:
        cmd = 'mkdir ' + folder
        subprocess.call(cmd, shell=True)
        return True


def replace_with_bases(base, read_result):
    """
    Parameters
    ----------
    base : string
        single character from the ref base ACGT.
    read_result : string
        line from the pileup read_result with commas and periods (,..).
    Returns
    -------
    replace commas and periods from results using the ref base
        DESCRIPTION.
    """
    if re.search("^[ACTG]", base):
        return read_result.replace(",", base[0]).replace(".", base[0])
    else:
        return read_result


def extract_haplogroups(path_markerfile, reads_thresh, base_majority,
                        path_pileupfile, log_output, fmf_output, outputfile, flag, use_old):
    print("Extracting haplogroups...")
    markerfile = pd.read_csv(path_markerfile, header=None, sep="\t")
    markerfile.columns = ["chr", "marker_name", "haplogroup", "pos", "mutation", "anc", "der"]
    markerfile = markerfile.drop_duplicates(subset='pos', keep='first', inplace=False)

    # packagemanagement is the best
    try:
        pileupfile = pd.read_csv(path_pileupfile, header=None, sep="\t",
                                 dtype={0: str, 1: int, 2: str, 3: int, 4: str, 5: str},
                                 on_bad_lines='skip')
    except TypeError:
        pileupfile = pd.read_csv(path_pileupfile, header=None, sep="\t",
                                 dtype={0: str, 1: int, 2: str, 3: int, 4: str, 5: str},
                                 error_bad_lines=False)

    pileupfile.columns = ['chr', 'pos', 'refbase', 'reads', 'align', 'quality']

    if flag == "cram":
        ref_base = pileupfile["refbase"].values
        read_results = pileupfile["align"].values
        new_read_results = list(map(replace_with_bases, ref_base, read_results))
        pileupfile["align"] = new_read_results

    log_output_list = ["Total of reads: " + str(len(pileupfile))]

    intersect_pos = np.intersect1d(pileupfile['pos'], markerfile['pos'])
    markerfile = markerfile.loc[markerfile['pos'].isin(intersect_pos)]
    markerfile = markerfile.sort_values(by=['pos'])
    pileupfile = pileupfile.loc[pileupfile['pos'].isin(intersect_pos)]

    pileupfile = pileupfile.drop(['chr'], axis=1)
    df = pd.merge(markerfile, pileupfile, on='pos')

    markerfile_len = len(markerfile)
    del [[pileupfile, markerfile]]
    gc.collect()

    # valid markers from positionsfile.txt
    log_output_list.append("Valid markers: " + str(markerfile_len))

    index_belowzero = df[df["reads"] == 0].index
    df_belowzero = df[df.index.isin(index_belowzero)]
    df_belowzero = df_belowzero.drop(['refbase', 'align', 'quality'], axis=1)
    df_belowzero["called_perc"] = "NA"
    df_belowzero["called_base"] = "NA"
    df_belowzero["state"] = "NA"
    df_belowzero["Description"] = "Position with zero reads"

    df = df[~df.index.isin(index_belowzero)]

    df_freq_table = get_frequency_table(df)
    df_freq_table = df_freq_table.drop(['+', '-'], axis=1)
    df = df.drop(['refbase', 'align', 'quality'], axis=1)

    list_col_indices = np.argmax(df_freq_table.values, axis=1)
    called_base = df_freq_table.columns[list_col_indices]
    total_count_bases = np.sum(df_freq_table.values, axis=1)
    max_count_bases = np.max(df_freq_table, axis=1)
    called_perc = round((max_count_bases / total_count_bases) * 100, 1)

    bool_anc = np.equal(np.array(called_base), df["anc"].values)
    bool_der = np.equal(np.array(called_base), df["der"].values)

    bool_list_anc = np.where(bool_anc, 'A', 'D')
    bool_list_anc = bool_list_anc.astype('object')
    bool_list_der = np.where(bool_der, 'D', 'A')
    bool_list_der = bool_list_der.astype('object')
    bool_list_state = np.equal(bool_list_anc, bool_list_der)

    df["called_perc"] = np.array(called_perc, dtype=int)
    df["called_base"] = called_base
    df["state"] = bool_list_anc
    df["bool_state"] = bool_list_state

    # discordant genotypes
    df_discordantgenotype = df[~bool_list_state]
    df_discordantgenotype = df_discordantgenotype.drop(["bool_state"], axis=1)
    df_discordantgenotype["state"] = "NA"
    df_discordantgenotype["Description"] = "Discordant genotype"
    columns_fmf = df_discordantgenotype.columns
    df = df[bool_list_state]

    # read threshold
    df_readsthreshold = df[df["reads"] < reads_thresh]
    df_readsthreshold["Description"] = "Below read threshold"
    df = df[df["reads"] >= reads_thresh]

    # filter by base percentage
    df_basemajority = df[df["called_perc"] < base_majority]
    df_basemajority["Description"] = "Below base majority"
    df = df[df["called_perc"] >= base_majority]

    df_fmf = pd.concat([df_belowzero, df_readsthreshold, df_basemajority, df_discordantgenotype], axis=0, sort=True)
    df_fmf = df_fmf[['chr', 'pos', 'marker_name', 'haplogroup', 'mutation', 'anc', 'der', 'reads',
                     'called_perc', 'called_base', 'state', 'Description']]

    df_out = df.drop(["bool_state"], axis=1)

    log_output_list.append("Markers with zero reads: " + str(len(df_belowzero)))
    log_output_list.append(
        "Markers below the read threshold {" + str(reads_thresh) + "}: " + str(len(df_readsthreshold)))
    log_output_list.append(
        "Markers below the base majority threshold {" + str(base_majority) + "}: " + str(len(df_basemajority)))
    log_output_list.append("Markers with discordant genotype: " + str(len(df_discordantgenotype)))
    log_output_list.append("Markers without haplogroup information: " + str(len(df_fmf)))
    log_output_list.append("Markers with haplogroup information: " + str(len(df_out)))

    with open(log_output, "a") as log:
        for marker in log_output_list:
            log.write(marker)
            log.write("\n")

    if use_old:
        df_out = df_out.sort_values(by=['haplogroup'], ascending=True)
        df_out = df_out[
            ["chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads", "called_perc", "called_base",
             "state"]]
        df_fmf.to_csv(fmf_output, sep="\t", index=False)
        df_out.to_csv(outputfile, sep="\t", index=False)
        return

    del [[df_basemajority, df_belowzero, df_discordantgenotype, df_readsthreshold, df_freq_table, df]]
    gc.collect()

    df_out = df_out[
        ["chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads", "called_perc", "called_base",
         "state"]]
    df_fmf.to_csv(fmf_output, sep="\t", index=False)

    # sort based on the tree
    lst_df = df_out.values.tolist()
    mappable_df = {}
    for lst in lst_df:
        if lst[3] not in mappable_df:
            mappable_df[lst[3]] = []
        mappable_df[lst[3]].append(lst)

    tree = Tree("Hg_Prediction_tables/tree.json")
    with open(outputfile, "w") as f:
        f.write('\t'.join(["chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads",
                           "called_perc", "called_base", "state", "depth\n"]))
        for node_key in tree.node_mapping:
            if node_key not in mappable_df:
                continue
            depth = tree.get(node_key).depth
            for lst in mappable_df[node_key]:
                f.write('\t'.join(map(str, lst)) + f"\t{depth}\n")


def samtools(threads, folder, folder_name, path_file, quality_thresh, markerfile, reference, flag, args, whole_time,
             use_old):
    file_name = folder_name
    outputfile = folder + "/" + folder_name + ".out"
    log_output = folder + "/" + folder_name + ".log"
    fmf_output = folder + "/" + folder_name + ".fmf"
    pileupfile = folder + "/" + folder_name + ".pu"

    start_time = time.time()
    if flag == "bam":
        if not os.path.exists(path_file + '.bai'):
            cmd = "samtools index -@{} {}".format(threads, path_file)
            print(cmd)
            subprocess.call(cmd, shell=True)
    elif flag == "cram":
        if not os.path.exists(path_file + '.crai'):
            cmd = "samtools index -@{} {}".format(threads, path_file)
            print(cmd)
            subprocess.call(cmd, shell=True)

    header, total_reads = chromosome_table(path_file, folder, file_name)

    execute_mpileup(header, path_file, pileupfile, quality_thresh, folder, reference)
    print("--- %.2f seconds in run PileUp ---" % (time.time() - start_time))

    start_time = time.time()
    extract_haplogroups(markerfile, args.Reads_thresh, args.Base_majority,
                        pileupfile, log_output, fmf_output, outputfile, flag, use_old)

    cmd = "rm {};".format(pileupfile)
    subprocess.call(cmd, shell=True)

    print("--- %.2f seconds in extracting haplogroups --- " % (time.time() - start_time))
    print("--- %.2f seconds to run Yleaf  ---" % (time.time() - whole_time))

    return outputfile


def logo():
    print(r"""

           |
          /|\          
         /\|/\    
        \\\|///   
         \\|//  
          |||   
          |||    
          |||    

        """)


def predict_haplogroup(source, path_file, output, use_old):
    if use_old:
        script = source + "/old_predict_haplogroup.py"
    else:
        script = source + "/predict_haplogroup.py"
    cmd = "python {} -i {} -o {}".format(script, path_file, output)
    subprocess.call(cmd, shell=True)


def main():
    whole_time = time.time()
    print(
        f"""\tErasmus MC Department of Genetic Identification \n\n\tYleaf: software tool for human Y-chromosomal
         \n\tphylogenetic analysis and haplogroup inference v{VERSION}\n""")
    logo()
    args = get_arguments()

    if args.use_old and "new" in args.position:
        print("ERROR: using the old prediction with new position files will give inaccurate predictions. Use files "
              "without new in the name.")
        sys.exit(-1)
    elif not args.use_old and "new" not in args.position:
        print("ERROR: using the new prediction method with old position files will give inaccurate predictions. Use"
              " files with new in the name.")
        sys.exit(-1)

    app_folder = os.path.dirname(os.path.realpath(__name__))
    out_path = args.Outputfile
    source = os.path.abspath(os.path.dirname(sys.argv[0]))
    out_folder = out_path
    hg_out = "hg_prediction.hg"
    if create_tmp_dirs(out_folder):
        if args.Fastq:
            files = check_if_folder(args.Fastq, '.fastq')
            for path_file in files:
                print(args.reference)
                if args.reference is None:
                    raise FileNotFoundError("-f missing reference")
                print("Starting...")
                folder_name = get_folder_name(path_file)
                folder = os.path.join(app_folder, out_folder, folder_name)
                if create_tmp_dirs(folder):
                    start_time = time.time()
                    sam_file = folder + "/" + folder_name + ".sam"
                    fastq_cmd = "bwa mem -t {} {} {} > {}".format(args.threads, args.reference, path_file, sam_file)
                    print(fastq_cmd)
                    subprocess.call(fastq_cmd, shell=True)
                    print("--- %s seconds in Indexing reads to reference ---" % (time.time() - start_time))
                    start_time = time.time()
                    bam_file = folder + "/" + folder_name + ".bam"
                    cmd = "samtools view -@ {} -bS {} | samtools sort -@ {} -m 2G -o {}".format(args.threads, sam_file,
                                                                                                args.threads, bam_file)
                    print(cmd)
                    subprocess.call(cmd, shell=True)
                    print("--- %s seconds in convertin Sam to Bam ---" % (time.time() - start_time))
                    cmd = "samtools index -@ {} {}".format(args.threads, bam_file)
                    subprocess.call(cmd, shell=True)
                    output_file = samtools(args.threads, folder, folder_name, bam_file, args.Quality_thresh,
                                           args.position, False, "bam", args, whole_time, args.use_old)
                    cmd = "rm {}".format(sam_file)
                    subprocess.call(cmd, shell=True)
            hg_out = out_folder + "/" + hg_out
            predict_haplogroup(source, out_folder, hg_out, args.use_old)
        elif args.Bamfile:
            files = check_if_folder(args.Bamfile, '.bam')
            for path_file in files:
                print("Starting...")
                print(path_file)
                bam_file = path_file
                folder_name = get_folder_name(path_file)
                folder = os.path.join(app_folder, out_folder, folder_name)
                if create_tmp_dirs(folder):
                    output_file = samtools(args.threads, folder, folder_name, bam_file, args.Quality_thresh,
                                           args.position, False, "bam", args, whole_time, args.use_old)
            hg_out = out_folder + "/" + hg_out
            predict_haplogroup(source, out_folder, hg_out, args.use_old)
        elif args.Cramfile:
            if args.reference is None:
                raise FileNotFoundError("-f missing reference")
            files = check_if_folder(args.Cramfile, '.cram')
            for path_file in files:
                print("Starting...")
                print(path_file)
                cram_file = path_file
                folder_name = get_folder_name(path_file)
                folder = os.path.join(app_folder, out_folder, folder_name)
                if create_tmp_dirs(folder):
                    output_file = samtools(args.threads, folder, folder_name, cram_file, args.Quality_thresh,
                                           args.position, args.reference, "cram", args, whole_time, args.use_old)
            hg_out = out_folder + "/" + hg_out
            predict_haplogroup(source, out_folder, hg_out, args.use_old)
    else:
        print("--- Yleaf failed! please check inputs... ---")


if __name__ == "__main__":
    main()
