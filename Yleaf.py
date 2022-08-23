#!/usr/bin/env python

"""
Yleaf detection of Y-Haplogroups in Human DNA v3.1

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Diego Montiel Gonzalez
Modified by: Bram van Wersch
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
from typing import Union, List, TextIO, Tuple, Dict
from collections import defaultdict

from tree import Tree

pd.options.mode.chained_assignment = None  # default='warn'

VERSION = 3.1

# this is not ideal
GENERAL_INFO = []

PREDICTION_OUT_FILE_NAME: str = "hg_prediction.hg"
HAPLOGROUP_IMAGE_FILE_NAME: str = "hg_tree_image.pdf"

LOG = logging.getLogger("yleaf_logger")


def main():

    print("Erasmus MC Department of Genetic Identification\nYleaf: software tool for human Y-chromosomal "
          f"phylogenetic analysis and haplogroup inference v{VERSION}")
    logo()

    args = get_arguments()
    out_folder = Path(args.output)
    safe_create_dir(out_folder, args.force)
    setup_logger(out_folder)

    LOG.info(f"Running Yleaf with command: {' '.join(sys.argv)}")

    app_folder = Path(__file__).absolute().parent
    source = app_folder
    if args.include_private_mutations and args.snp_file is None:
        LOG.error("Missing snp_file. Please specify one when requesting to include private mutations. These files can"
                  " be found in the snp_data_tables folder.")
        raise SystemExit("Missing snp_file argument. Please supply it.")
    if args.fastq:
        main_fastq(args, app_folder, out_folder)
    elif args.bamfile:
        main_bam_cram(args, app_folder, out_folder, True)
    elif args.cramfile:
        main_bam_cram(args, app_folder, out_folder, False)
    else:
        LOG.error("Please specify either a bam, a cram or a fastq file")
        raise ValueError("Please specify either a bam, a cram or a fastq file")
    hg_out = out_folder / PREDICTION_OUT_FILE_NAME
    predict_haplogroup(source, out_folder, hg_out, args.use_old)
    if args.draw_haplogroups:
        draw_haplogroups(hg_out, args.collapsed_draw_mode)

    LOG.info("Done!")


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
                        help="fasta reference genome sequence ", metavar="PATH", required=False, type=check_file,
                        default=None)

    parser.add_argument("-pos", "--position",
                        help="Positions file found in the Position_files folder. Use old_position_files when using the "
                             " -old flag, otherwise use the new_position_files.", metavar="PATH", required=True,
                        type=check_file)

    parser.add_argument("-out", "--output", required=True,
                        help="Folder name containing outputs", metavar="STRING")

    parser.add_argument("-r", "--reads_treshold",
                        help="The minimum number of reads for each base. (default=50)",
                        type=int, required=False,
                        default=50, metavar="INT")

    parser.add_argument("-q", "--quality_thresh",
                        help="Minimum quality for each read, integer between 10 and 40. [10-40]",
                        type=int, required=True, metavar="INT")

    parser.add_argument("-b", "--base_majority",
                        help="The minimum percentage of a base result for acceptance, integer between 50 and 99."
                             " [50-99]",
                        type=int, required=True, metavar="INT")

    parser.add_argument("-t", "--threads", dest="threads",
                        help="Set number of additional threads to use during alignment BWA-MEM",
                        type=int, default=1, metavar="INT")

    parser.add_argument("-old", "--use_old", dest="use_old",
                        help="Add this value if you want to use the old prediction method of Yleaf (version 2.3). This"
                             " version only uses the ISOGG tree and slightly different prediction criteria.",
                        action="store_true")
    parser.add_argument("-dh", "--draw_haplogroups", help="Draw the predicted haplogroups in the haplogroup tree.",
                        action="store_true")
    parser.add_argument("-hc", "--collapsed_draw_mode",
                        help="Add this flag to compress the haplogroup tree image and remove all uninformative "
                             "haplogroups from it.", action="store_true")
    # TODO improve help
    parser.add_argument("-pm", "--include_private_mutations", help="Include private mutations. These are mutations not"
                                                                   " captured in a haplogroup but observed in the "
                                                                   "individual.",
                        action="store_true")
    parser.add_argument("-sf", "--snp_file", help="File containing snp's with confirmed rates, for the identification"
                                                  " of private mutations",
                        metavar="FILE", type=check_file, default=None)
    parser.add_argument("-mr", "--maximum_mutation_rate", help="Maximum rate of minor allele for it to be considered"
                                                               " as a private mutation. (default=0.01)",
                        default=0.01, type=float, metavar="FLOAT")

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


def setup_logger(
    out_folder: Path
):
    """Setup logging"""
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(levelname)s (%(relativeCreated)d ms) - %(message)s')
    handler.setFormatter(formatter)
    LOG.addHandler(handler)

    file_handler = logging.FileHandler(filename=out_folder / "run.log")
    file_handler.setFormatter(formatter)
    LOG.addHandler(file_handler)

    LOG.setLevel(logging.INFO)
    LOG.debug("Logger created")


def safe_create_dir(
    folder: Path,
    force: bool
):
    """Create the given folder. If the folder is already present delete if the user agrees."""
    if folder.is_dir():
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
            print("Failed to create directory. Exiting...")
            raise


def call_command(
        command_str: str,
        stdout_location: TextIO = None
):
    """Call a command on the command line and make sure to exit if it fails."""
    LOG.info(f"Started running the following command: {command_str}")
    if stdout_location is None:
        process = subprocess.Popen(command_str, stderr=subprocess.PIPE, shell=True)
    else:
        process = subprocess.Popen(command_str, stderr=subprocess.PIPE, stdout=stdout_location, shell=True)
    # blocking call
    stdout, stderr = process.communicate()
    # will only fail if returncode is not 0
    if process.returncode != 0:
        LOG.error(f"Above call failed with message {stderr.decode('utf-8')}")
        raise SystemExit("Failed command execution")
    LOG.info("Finished running the command")


def main_fastq(
    args: argparse.Namespace,
    app_folder: Path,
    out_folder: Path
):
    if args.reference is None:
        LOG.error("Missing reference genome, please supply one with -f")
        raise SystemExit("Missing reference genome, please supply one with -f")
    files = get_files_with_extension(args.fastq, '.fastq')
    for path_file in files:
        LOG.info(f"Starting with running for {path_file}")
        folder = app_folder / out_folder / path_file.name.split(".")[0]
        safe_create_dir(folder, args.force)
        sam_file = folder / (folder.name + ".sam")

        fastq_cmd = f"minimap2 -ax sr -t {args.threads} {args.reference} {path_file} > {sam_file}"
        call_command(fastq_cmd)
        bam_file = folder / (folder.name + ".bam")
        cmd = "samtools view -@ {} -bS {} | samtools sort -@ {} -m 2G -o {}".format(args.threads, sam_file,
                                                                                    args.threads, bam_file)
        call_command(cmd)
        cmd = "samtools index -@ {} {}".format(args.threads, bam_file)
        call_command(cmd)
        samtools(folder, bam_file, None, True, args)
        os.remove(sam_file)
        if args.include_private_mutations:
            find_private_mutations(folder, path_file, reference, args)
        write_info_file(folder, folder.name)
        LOG.info(f"Finished running for {path_file.name}")
        print()


def main_bam_cram(
    args: argparse.Namespace,
    app_folder: Path,
    out_folder: Path,
    is_bam: bool
):
    if args.reference is None and not is_bam:
        LOG.error("Missing reference genome, please supply one with -f")
        raise SystemExit("Missing reference genome, please supply one with -f")
    if is_bam:
        reference = None
    else:
        reference = args.reference
    files = get_files_with_extension(args.bamfile, '.bam')
    LOG.info("Starting with running for bamfile...")
    for path_file in files:
        LOG.info(f"Starting with running for {path_file}")
        folder = app_folder / out_folder / path_file.name.split(".")[0]
        safe_create_dir(folder, args.force)
        samtools(folder, path_file, reference, is_bam, args)
        if args.include_private_mutations:
            find_private_mutations(folder, path_file, reference, args)
        write_info_file(folder, folder.name)
        LOG.info(f"Finished running for {path_file.name}")
        print()


def find_private_mutations(
    output_folder: Path,
    path_file: Path,
    reference: Union[Path, None],
    args: argparse.Namespace
):
    LOG.info("Starting with extracting private mutations...")

    # run mpileup
    pileup_file = output_folder / "temp_pileup.pu"
    execute_mpileup(None, path_file, pileup_file, args.quality_thresh, reference)

    # get positions with known SNPs that are not part of a haplogroup
    filter_positions = set()
    with open(args.position) as f:
        for line in f:
            filter_positions.add(line.strip().split("\t")[3])

    positions_of_interest = defaultdict(list)
    with open(args.snp_file) as f:
        f.readline()
        for line in f:
            rs, position, major_allele, minor_allele, frequency = line.strip().split(",")
            if float(frequency) > args.maximum_mutation_rate:
                continue
            positions_of_interest[position].append({"rs": rs, "major_allele": major_allele,
                                                    "minor_allele": minor_allele, "frequency": frequency})
    positions_of_interest = dict(positions_of_interest)
    private_mutations = []
    with open(pileup_file) as f:
        for index, line in enumerate(f):
            try:
                chrom, position, ref_base, count, aligned, quality = line.strip().split()
            except ValueError:
                LOG.warning(f"Failed to read line {index} of pileupfile")
                continue
            # not on y-chromosome
            if "y" not in chrom.split("_")[0].lower():
                continue
            # not enough reads
            if int(count) < args.reads_treshold:
                continue
            if reference is not None:
                aligned = replace_with_bases(ref_base, aligned)
            if position in filter_positions:
                continue
            freq_dict = get_frequencies(aligned)
            actual_allele, allele_count = max(freq_dict.items(), key=lambda x: x[1])
            called_percentage = allele_count / sum(freq_dict.values()) * 100
            # not a high enough base majority measured, cannot be sure of real allele
            if called_percentage < args.base_majority:
                continue
            # not annotated in dbsnp
            if position not in positions_of_interest:
                continue
                # private_mutations.append(f"{chrom}\t{position}\t-\t{ref_base}->{actual_allele}\t{ref_base}\t"
                #                          f"{actual_allele}\t{allele_count}\t{called_percentage}\t{frequency}\n")
            else:
                possible_minor_alleles = {dct["minor_allele"] for dct in positions_of_interest[position]}
                # measured allele is not the dbsnp allele
                if actual_allele not in possible_minor_alleles:
                    continue

                matched_pos_dct = [dct for dct in positions_of_interest[position]
                                   if dct["minor_allele"] == actual_allele][0]
                rs, major_allele, minor_allele, frequency = matched_pos_dct.values()
                private_mutations.append(f"{chrom}\t{position}\t{rs}\t{major_allele}->{actual_allele}\t{major_allele}\t"
                                         f"{actual_allele}\t{allele_count}\t{called_percentage}\t{frequency}\n")
    os.remove(pileup_file)
    with open(output_folder / "private_mutations.csv", "w") as f:
        f.write(''.join(private_mutations))
    LOG.info("Finished extracting private mutations")


def get_frequency_table(
    mpileup: List[str]
) -> Dict[str, List[int]]:
    frequency_table = {}
    for i in mpileup:
        fastadict = get_frequencies(i[9])
        frequency_table[i[3]] = list(fastadict.values())
    return frequency_table


def get_frequencies(
    sequence: str
) -> Dict[str, int]:
    fastadict = {"A": 0, "T": 0, "G": 0, "C": 0}
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
    return fastadict


def find_all_indels(
    s: str
) -> List[int]:
    list_pos = []
    for i in find_all(s, "-"):
        list_pos.append(i)
    for i in find_all(s, "+"):
        list_pos.append(i)
    return sorted(list_pos)


def count_indels(
    s: str,
    pos: List[int]
) -> Dict[str, int]:
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


def trimm_indels(
    s: str,
    pos: List[int]
) -> str:
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


def trimm_caret(
    s: str
) -> str:
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


def find_all(
    c: str,
    s: str
) -> List[int]:
    return [x for x in range(c.find(s), len(c)) if c[x] == s]


def write_bed_file(
    bed: Path,
    markerfile: Path,
    header: str
):
    mf = pd.read_csv(markerfile, sep="\t", header=None)
    mf = mf[[0, 3]]
    mf[0] = header
    mf.to_csv(str(bed), sep="\t", index=False, header=False)


def execute_mpileup(
    bed: Union[Path, None],
    bam_file: Path,
    pileupfile: Path,
    quality_thresh: float,
    reference: Union[Path, None]
):
    cmd = "samtools mpileup"
    if bed is not None:
        cmd += f" -l {str(bed)}"

    if reference is not None:
        cmd += f" -f {str(reference)}"
    cmd += f" -AQ{quality_thresh} {str(bam_file)} > {str(pileupfile)}"
    call_command(cmd)


def chromosome_table(
    path_file: Path,
    path_folder: Path,
    file_name: str
) -> Tuple[str, int, int]:
    output = path_folder / (file_name + '.chr')
    tmp_output = path_folder / "tmp.txt"
    f = open(tmp_output, "w")
    cmd = f"samtools idxstats {path_file}"
    call_command(cmd, stdout_location=f)
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
    else:
        LOG.error("Unable to find Y-chromosomal data in the provided files. Exiting...")
        raise SystemExit("No Y-chromosomal data")


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


def replace_with_bases(
    base: str,
    read_result: str
) -> str:
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
    use_old: bool
):
    LOG.info("Starting with extracting haplogroups...")
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

    intersect_pos = np.intersect1d(pileupfile['pos'], markerfile['pos'])
    markerfile = markerfile.loc[markerfile['pos'].isin(intersect_pos)]
    markerfile = markerfile.sort_values(by=['pos'])
    pileupfile = pileupfile.loc[pileupfile['pos'].isin(intersect_pos)]

    pileupfile = pileupfile.drop(['chr'], axis=1)
    df = pd.merge(markerfile, pileupfile, on='pos')

    markerfile_len = len(markerfile)

    # valid markers from positionsfile.txt
    GENERAL_INFO.append("Valid markers: " + str(markerfile_len))

    index_belowzero = df[df["reads"] == 0].index
    df_belowzero = df[df.index.isin(index_belowzero)]
    df_belowzero = df_belowzero.drop(['refbase', 'align', 'quality'], axis=1)
    df_belowzero["called_perc"] = "NA"
    df_belowzero["called_base"] = "NA"
    df_belowzero["state"] = "NA"
    df_belowzero["Description"] = "Position with zero reads"

    df = df[~df.index.isin(index_belowzero)]

    freq_dict = get_frequency_table(df.values)
    df_freq_table = pd.DataFrame.from_dict(freq_dict, orient='index')
    df_freq_table.columns = ["A", "T", "G", "C", "+", "-"]

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

    GENERAL_INFO.append("Markers with zero reads: " + str(len(df_belowzero)))
    GENERAL_INFO.append(
        "Markers below the read threshold {" + str(reads_thresh) + "}: " + str(len(df_readsthreshold)))
    GENERAL_INFO.append(
        "Markers below the base majority threshold {" + str(base_majority) + "}: " + str(len(df_basemajority)))
    GENERAL_INFO.append("Markers with discordant genotype: " + str(len(df_discordantgenotype)))
    GENERAL_INFO.append("Markers without haplogroup information: " + str(len(df_fmf)))
    GENERAL_INFO.append("Markers with haplogroup information: " + str(len(df_out)))

    if use_old:
        df_out = df_out.sort_values(by=['haplogroup'], ascending=True)
        df_out = df_out[
            ["chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads", "called_perc", "called_base",
             "state"]]
        df_fmf.to_csv(fmf_output, sep="\t", index=False)
        df_out.to_csv(outputfile, sep="\t", index=False)
        return

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
    output_folder: Path,
    path_file: Path,
    reference: Union[Path, None],
    is_bam_pathfile: bool,
    args: argparse.Namespace,
):

    outputfile = output_folder / (output_folder.name + ".out")
    fmf_output = output_folder / (output_folder.name + ".fmf")
    pileupfile = output_folder / "temp_pileup.pu"

    if is_bam_pathfile:
        if not any([Path(str(path_file) + ".bai").exists(), Path(str(path_file).rstrip(".bam") + '.bai').exists()]):
            cmd = "samtools index -@{} {}".format(args.threads, path_file)
            call_command(cmd)
    else:
        if not any([Path(str(path_file) + ".crai").exists(), Path(str(path_file).rstrip(".cram") + '.crai').exists()]):
            cmd = "samtools index -@{} {}".format(args.threads, path_file)
            call_command(cmd)
    header, mapped, unmapped = chromosome_table(path_file, output_folder, output_folder.name)

    bed = Path(str(args.position).rsplit(".", 1)[0] + ".bed")
    if not bed.exists():
        write_bed_file(bed, args.position, header)

    execute_mpileup(bed, path_file, pileupfile, args.quality_thresh, reference)

    GENERAL_INFO.append("Total of mapped reads: " + str(mapped))
    GENERAL_INFO.append("Total of unmapped reads: " + str(unmapped))

    extract_haplogroups(args.position, args.reads_treshold, args.base_majority,
                        pileupfile, fmf_output, outputfile, is_bam_pathfile, args.use_old)

    os.remove(pileupfile)

    LOG.info("Finished extracting haplogroups")


def write_info_file(
    folder: Path,
    folder_name: str
):
    # Write all information in GENERAL_INFO
    global GENERAL_INFO
    try:
        with open(folder / (folder_name + ".info"), "a") as f:
            for marker in GENERAL_INFO:
                f.write(marker)
                f.write("\n")
    except IOError:
        LOG.warning("Failed to write .info file")
    GENERAL_INFO = [GENERAL_INFO[0]]  # make sure to reset except for the command used to call


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
        cmd = "python {} -i {} -o {}".format(script, path_file, output)
        call_command(cmd)
    else:
        import predict_haplogroup
        namespace = argparse.Namespace(input=path_file, outfile=output,
                                       minimum_score=predict_haplogroup.DEFAULT_MIN_SCORE)
        predict_haplogroup.main(namespace)


def draw_haplogroups(
    haplogroup_file: Path,
    collapsed_draw_mode: bool
):
    # make sure that it is only imported if requested by user
    import haplogroup_tree_image
    namespace = argparse.Namespace(input=haplogroup_file, collapse_mode=collapsed_draw_mode,
                                   outfile=haplogroup_file.parent / HAPLOGROUP_IMAGE_FILE_NAME)
    haplogroup_tree_image.main(namespace)


if __name__ == "__main__":
    main()
