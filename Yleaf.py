#!/usr/bin/env python

"""
Yleaf detection of Y-Haplogroups in Human DNA v3.0

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Diego Montiel Gonzalez
Slightly modified by: Bram van Wersch
"""
import argparse
import os
import sys
import re
import time
import logging
import shutil
import subprocess
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from pathlib import Path
from typing import Union, List

from tree import Tree

pd.options.mode.chained_assignment = None  # default='warn'

VERSION = 3.0

# this is not ideal
LOG_LIST = []

PREDICTION_OUT_FILE_NAME: str = "hg_prediction.hg"

LOG = logging.getLogger("yleaf_logger")


def main():

    setup_logger()
    LOG.info("Erasmus MC Department of Genetic Identification\nYleaf: software tool for human Y-chromosomal "
             f"phylogenetic analysis and haplogroup inference v{VERSION}")
    logo()
    args = get_arguments()

    LOG_LIST.append("Command: " + ' '.join(sys.argv))

    app_folder = Path(__file__).absolute().parent
    out_folder = Path(args.output)
    source = app_folder
    safe_create_dir(out_folder, args.force)
    LOG.info("Starting with Yleaf run...")
    if args.fastq:
        main_fastq(args, app_folder, out_folder, source)
    elif args.bamfile:
        main_bam(args, app_folder, out_folder, source)
    elif args.cramfile:
        if args.reference is None:
            raise FileNotFoundError("-f missing reference")
        files = get_files_with_extension(args.cramfile, '.cram')
        for path_file in files:
            print(path_file)
            cram_file = path_file
            folder_name = get_folder_name(path_file)
            folder = app_folder / out_folder / folder_name
            safe_create_dir(folder, args.force)
            samtools(args.threads, folder, folder_name, cram_file, args.Quality_thresh, args.position,
                     args.reference, False, args, args.use_old)
            write_log_file(folder, folder_name)
        hg_out = out_folder / PREDICTION_OUT_FILE_NAME
        predict_haplogroup(source, out_folder, hg_out, args.use_old)
    else:
        raise ValueError("Please specify either a bam, a cram or a fastq file")
    LOG.info("Done!")


def setup_logger():
    """Setup logging"""
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(levelname)s (%(relativeCreated)d ms) - %(message)s')
    handler.setFormatter(formatter)
    LOG.addHandler(handler)
    LOG.setLevel(logging.DEBUG)
    LOG.debug("Logger created")


def get_arguments() -> argparse.Namespace:
    parser = ArgumentParser()

    parser.add_argument("-fastq", "--fastq", required=False,
                        help="Use raw FastQ files", metavar="PATH", type=check_file)

    parser.add_argument("-bam", "--bamfile", required=False,
                        help="input BAM file", metavar="PATH", type=check_file)

    parser.add_argument("-cram", "--cramfile", required=False,
                        help="input CRAM file", metavar="PATH", type=check_file)

    parser.add_argument("-force", "--force", action="store_true",
                        help="Delete files without asking")

    parser.add_argument("-f", "--reference",
                        help="fasta reference genome sequence ", metavar="PATH", required=False, type=check_file)

    parser.add_argument("-pos", "--position",
                        help="Positions file found in the Position_files folder. Use old_position_files when using the "
                             " -old flag, otherwise use the new_position_files.", metavar="PATH", required=True,
                        type=check_file)

    parser.add_argument("-out", "--output", required=True,
                        help="Folder name containing outputs", metavar="STRING")

    parser.add_argument("-r", "--reads_treshold",
                        help="The minimum number of reads for each base. (default=50)",
                        type=int, required=False,
                        default=50)

    parser.add_argument("-q", "--Quality_thresh",
                        help="Minimum quality for each read, integer between 10 and 40. [10-40]",
                        type=int, required=True)

    parser.add_argument("-b", "--Base_majority",
                        help="The minimum percentage of a base result for acceptance, integer between 50 and 99."
                             " [50-99]",
                        type=int, required=True)

    parser.add_argument("-t", "--Threads", dest="threads",
                        help="Set number of additional threads to use during alignment BWA-MEM",
                        type=int,
                        default=1)

    parser.add_argument("-old", "--use_old", dest="use_old",
                        help="Add this value if you want to use the old prediction method of Yleaf (version 2.3). This"
                             " version only uses the ISOGG tree and slightly different prediction criteria.",
                        action="store_true")

    args = parser.parse_args()
    return args


def check_file(
    path: str
) -> Path:
    """Check for the presence of a file and return a Path object"""
    object_path = Path(path)
    if not object_path.exists():
        raise argparse.ArgumentTypeError("Path to provided file/dir does not exist")
    return object_path


def safe_create_dir(
    folder: Union[str, Path],
    force: bool
):
    """Create the given folder. If the folder is already present delete if the user agrees."""
    if os.path.isdir(folder):
        while True and not force:
            LOG.warning("Folder " + str(folder) + " already exists, would you like to remove it?")
            choice = input("y/n: ")
            if str(choice).upper() == "Y":
                break
            elif str(choice).upper() == "N":
                sys.exit(0)
            else:
                print("Please type y/Y or n/N")
        shutil.rmtree(folder)
        os.mkdir(folder)
    else:
        try:
            os.mkdir(folder)
        except OSError:
            LOG.error("WARNING: failed to create directory. Exiting...")
            raise


def main_fastq(
    args: argparse.Namespace,
    app_folder: Path,
    out_folder: Path,
    source: Path
):
    if args.reference is None:
        LOG.error("Missing reference genome, please supply one with -f")
        raise SystemExit("Missing reference genome, please supply one with -f")
    files = get_files_with_extension(args.fastq, '.fastq')
    for path_file in files:
        folder_name = get_folder_name(path_file)
        folder = app_folder / out_folder / folder_name
        safe_create_dir(folder, args.force)
        sam_file = folder / (folder_name + ".sam")
        fastq_cmd = "./bwa-mem2 mem -t {} {} {} > {}".format(args.threads, args.reference, path_file, sam_file)
        call_command(fastq_cmd)
        bam_file = folder / (folder_name + ".bam")
        cmd = "samtools view -@ {} -bS {} | samtools sort -@ {} -m 2G -o {}".format(args.threads, sam_file,
                                                                                    args.threads, bam_file)
        call_command(cmd)
        cmd = "samtools index -@ {} {}".format(args.threads, bam_file)
        call_command(cmd)
        samtools(args.threads, folder, folder_name, bam_file, args.Quality_thresh, args.position, False,
                 True, args, args.use_old)
        os.remove(sam_file)
        write_log_file(folder, folder_name)
    hg_out = out_folder / PREDICTION_OUT_FILE_NAME
    predict_haplogroup(source, out_folder, hg_out, args.use_old)


def main_bam(args, app_folder, out_folder, source):
    files = get_files_with_extension(args.bamfile, '.bam')
    LOG.info("Starting with running for bamfile...")
    for path_file in files:
        bam_file = path_file
        folder_name = get_folder_name(path_file)
        folder = app_folder / out_folder / folder_name
        safe_create_dir(folder, args.force)
        samtools(args.threads, folder, folder_name, bam_file, args.Quality_thresh, args.position, False,
                 True, args, args.use_old)
        write_log_file(folder, folder_name)
        LOG.info(f"Finished running for {path_file.name}")
    hg_out = out_folder / PREDICTION_OUT_FILE_NAME
    predict_haplogroup(source, out_folder, hg_out, args.use_old)


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


def call_command(command_str: str):
    LOG.info(f"Started running the following command: {command_str}")
    process = subprocess.Popen(command_str, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    # blocking call
    stdout, stderr = process.communicate()
    # will only fail if returncode is not 0
    if process.returncode != 0:
        LOG.error(f"Above call failed with message {stderr.decode('utf-8')}")
        raise SystemExit("Failed command execution")
    LOG.info("Finished running the command")


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


def check_bed(bed, markerfile, header):
    if not os.path.isfile(bed):
        mf = pd.read_csv(markerfile, sep="\t", header=None)
        mf = mf[[0, 3]]
        mf[0] = header

        mf.to_csv(bed, sep="\t", index=False, header=False)


def execute_mpileup(bed, bam_file, pileupfile, quality_thresh, reference):

    if reference:
        cmd = "samtools mpileup -l {} -f {} -AQ{} {} > {}".format(bed, reference, quality_thresh, bam_file, pileupfile)
    else:
        cmd = "samtools mpileup -l {} -AQ{} {} > {}".format(bed, quality_thresh, bam_file, pileupfile)
    subprocess.call(cmd, shell=True)


def chromosome_table(path_file: Path, path_folder: Path, file_name: str):
    output = path_folder / (file_name + '.chr')
    tmp_output = path_folder / "tmp.txt"
    print()
    f = open(tmp_output, "w")
    subprocess.call(["samtools", "idxstats", str(path_file)], stdout=f)
    df_chromosome = pd.read_table(tmp_output, header=None)

    total_reads = sum(df_chromosome[2])

    unmapped = df_chromosome[df_chromosome[0].str.contains("Y")][3].values[0]

    df_chromosome["perc"] = (df_chromosome[2] / total_reads) * 100
    df_chromosome = df_chromosome.round(decimals=2)
    df_chromosome['perc'] = df_chromosome['perc'].astype(str) + '%'
    df_chromosome = df_chromosome.drop(columns=[1, 3])
    df_chromosome.columns = ['chr', 'reads', 'perc']
    df_chromosome.to_csv(output, index=None, sep="\t")

    f.close()
    os.remove(tmp_output)

    if 'Y' in df_chromosome["chr"].values:
        return "Y", total_reads, unmapped
    elif 'chrY' in df_chromosome["chr"].values:
        return "chrY", total_reads, unmapped


def get_files_with_extension(
    path: Union[str, Path],
    ext: str
) -> List[Path]:
    """Get all files with a certain extension from a path. The path can be a file or a dir."""
    filtered_files = []
    path = Path(path)  # to be sure
    if path.is_dir():
        for file in path.iterdir():
            if file.suffix == ext:
                filtered_files.append(file)
        return filtered_files
    else:
        return [path]


def get_folder_name(
    path_file: Path
) -> str:
    file = path_file.name
    folder_name = os.path.splitext(file)[0]
    return folder_name


def replace_with_bases(base, read_result):
    print(base, read_result)
    if re.search("^[ACTG]", base):
        return read_result.replace(",", base[0]).replace(".", base[0])
    else:
        return read_result


def extract_haplogroups(
    path_markerfile: Path,
    reads_thresh: float,
    base_majority: int,
    path_pileupfile: Path,
    fmf_output: Path,
    outputfile: Path,
    is_bam_file: bool,
    use_old: bool,
    mapped: int,
    unmapped: int
):
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

    if not is_bam_file:
        ref_base = pileupfile["refbase"].values
        read_results = pileupfile["align"].values
        new_read_results = list(map(replace_with_bases, ref_base, read_results))
        pileupfile["align"] = new_read_results

    LOG_LIST.append("Total of mapped reads: " + str(mapped))
    LOG_LIST.append("Total of unmapped reads: " + str(unmapped))

    intersect_pos = np.intersect1d(pileupfile['pos'], markerfile['pos'])
    markerfile = markerfile.loc[markerfile['pos'].isin(intersect_pos)]
    markerfile = markerfile.sort_values(by=['pos'])
    pileupfile = pileupfile.loc[pileupfile['pos'].isin(intersect_pos)]

    pileupfile = pileupfile.drop(['chr'], axis=1)
    df = pd.merge(markerfile, pileupfile, on='pos')

    markerfile_len = len(markerfile)
    # del [[pileupfile, markerfile]]
    # gc.collect()

    # valid markers from positionsfile.txt
    LOG_LIST.append("Valid markers: " + str(markerfile_len))

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

    LOG_LIST.append("Markers with zero reads: " + str(len(df_belowzero)))
    LOG_LIST.append(
        "Markers below the read threshold {" + str(reads_thresh) + "}: " + str(len(df_readsthreshold)))
    LOG_LIST.append(
        "Markers below the base majority threshold {" + str(base_majority) + "}: " + str(len(df_basemajority)))
    LOG_LIST.append("Markers with discordant genotype: " + str(len(df_discordantgenotype)))
    LOG_LIST.append("Markers without haplogroup information: " + str(len(df_fmf)))
    LOG_LIST.append("Markers with haplogroup information: " + str(len(df_out)))

    if use_old:
        df_out = df_out.sort_values(by=['haplogroup'], ascending=True)
        df_out = df_out[
            ["chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads", "called_perc", "called_base",
             "state"]]
        df_fmf.to_csv(fmf_output, sep="\t", index=False)
        df_out.to_csv(outputfile, sep="\t", index=False)
        return

    # del [[df_basemajority, df_belowzero, df_discordantgenotype, df_readsthreshold, df_freq_table, df]]
    # gc.collect()

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


def samtools(
    threads: int,
    folder: Path,
    folder_name: str,
    path_file: Path,
    quality_thresh: float,
    markerfile: Path,
    reference: bool,
    is_bam_pathfile: bool,
    args: argparse.Namespace,
    use_old: bool
):
    file_name = folder_name
    outputfile = folder / (folder_name + ".out")
    fmf_output = folder / (folder_name + ".fmf")
    pileupfile = folder / (folder_name + ".pu")

    if is_bam_pathfile:
        if not os.path.exists(str(path_file) + '.bai'):
            cmd = "samtools index -@{} {}".format(threads, path_file)
            LOG.info(f"Calling samtools with: {cmd}")
            subprocess.call(cmd, shell=True)
    else:
        if not os.path.exists(path_file + '.crai'):
            cmd = "samtools index -@{} {}".format(threads, path_file)
            LOG.info(f"Calling samtools with: {cmd}")
            subprocess.call(cmd, shell=True)
    header, mapped, unmapped = chromosome_table(path_file, folder, file_name)

    bed = str(markerfile).rsplit(".", 1)[0] + ".bed"
    check_bed(bed, markerfile, header)

    execute_mpileup(bed, path_file, pileupfile, quality_thresh, reference)
    LOG.info("Finished running PileUp")

    start_time = time.time()
    extract_haplogroups(markerfile, args.reads_treshold, args.Base_majority,
                        pileupfile, fmf_output, outputfile, is_bam_pathfile, use_old, mapped, unmapped)

    os.remove(pileupfile)
    os.remove(bed)

    print("--- %.2f seconds in extracting haplogroups --- " % (time.time() - start_time))


def write_log_file(
    folder: Path,
    folder_name: str
):
    # Write all information in LOG_LIST
    global LOG_LIST
    try:
        with open(folder / (folder_name + ".log"), "a") as log:
            for marker in LOG_LIST:
                log.write(marker)
                log.write("\n")
    except IOError:
        print("Failed to write .log file")
    LOG_LIST = [LOG_LIST[0]]  # make sure to reset except for the command used to call


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


def predict_haplogroup(
    source: Path,
    path_file: Path,
    output: Path,
    use_old: bool
):
    if use_old:
        script = source / "old_predict_haplogroup.py"
    else:
        script = source / "predict_haplogroup.py"
    cmd = "python {} -i {} -o {}".format(script, path_file, output)
    subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    main()
