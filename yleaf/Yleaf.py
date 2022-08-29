#!/usr/bin/env python

"""
Yleaf detection of Y-Haplogroups in Human DNA v3.1

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Diego Montiel Gonzalez
Extensively modified by: Bram van Wersch
"""

import argparse
import os
import sys
import logging
import shutil
import subprocess
import pandas as pd
import numpy as np
import multiprocessing
from functools import partial
from argparse import ArgumentParser
from pathlib import Path
from typing import Union, List, TextIO, Tuple, Dict, Set
from collections import defaultdict

from yleaf import __version__
from yleaf.tree import Tree
from yleaf import yleaf_constants, download_reference

pd.options.mode.chained_assignment = None  # default='warn'

CACHED_POSITION_DATA: Union[Set[str], None] = None
CACHED_SNP_DATABASE: Union[Dict[str, List[Dict[str, str]]], None] = None
CACHED_REFERENCE_FILE: Union[List[str], None] = None
NUM_SET: Set[str] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"}
ACCEPTED_REF_BASES: Set[str] = {"A", "C", "G", "T"}

# path constants
PREDICTION_OUT_FILE_NAME: str = "hg_prediction.hg"
HAPLOGROUP_IMAGE_FILE_NAME: str = "hg_tree_image"

LOG: logging.Logger = logging.getLogger("yleaf_logger")


def main():

    print("Erasmus MC Department of Genetic Identification\nYleaf: software tool for human Y-chromosomal "
          f"phylogenetic analysis and haplogroup inference v{__version__}")
    logo()

    args = get_arguments()
    out_folder = Path(args.output)
    safe_create_dir(out_folder, args.force)
    setup_logger(out_folder)

    LOG.info(f"Running Yleaf with command: {' '.join(sys.argv)}")

    if args.fastq:
        main_fastq(args, out_folder)
    elif args.bamfile:
        main_bam_cram(args, out_folder, True)
    elif args.cramfile:
        main_bam_cram(args, out_folder, False)
    else:
        LOG.error("Please specify either a bam, a cram or a fastq file")
        raise ValueError("Please specify either a bam, a cram or a fastq file")
    hg_out = out_folder / PREDICTION_OUT_FILE_NAME
    predict_haplogroup(out_folder, hg_out, args.use_old)
    if args.draw_haplogroups:
        draw_haplogroups(hg_out, args.collapsed_draw_mode)

    LOG.info("Done!")


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
    parser.add_argument("-rg", "--reference_genome",
                        help="The reference genome of your bam/cram or fastq file.",
                        choices=[yleaf_constants.HG19, yleaf_constants.HG38], required=True)
    parser.add_argument("-o", "--output", required=True,
                        help="Folder name containing outputs", metavar="STRING")
    parser.add_argument("-r", "--reads_treshold",
                        help="The minimum number of reads for each base. (default=50)",
                        type=int, required=False,
                        default=50, metavar="INT")
    parser.add_argument("-q", "--quality_thresh",
                        help="Minimum quality for each read, integer between 10 and 40. [10-40] (default=10)",
                        type=int, required=False, metavar="INT", default=10)
    parser.add_argument("-b", "--base_majority",
                        help="The minimum percentage of a base result for acceptance, integer between 50 and 99."
                             " [50-99] (default=50)",
                        type=int, required=False, metavar="INT", default=50)
    parser.add_argument("-t", "--threads", dest="threads",
                        help="The number of processes to use when running Yleaf.",
                        type=int, default=1, metavar="INT")

    # arguments for prediction
    parser.add_argument("-old", "--use_old", dest="use_old",
                        help="Add this value if you want to use the old prediction method of Yleaf (version 2.3). This"
                             " version only uses the ISOGG tree and slightly different prediction criteria.",
                        action="store_true")

    # arguments for drawing haplo group trees
    parser.add_argument("-dh", "--draw_haplogroups", help="Draw the predicted haplogroups in the haplogroup tree.",
                        action="store_true")
    parser.add_argument("-hc", "--collapsed_draw_mode",
                        help="Add this flag to compress the haplogroup tree image and remove all uninformative "
                             "haplogroups from it.", action="store_true")

    # arguments for private mutations
    parser.add_argument("-p", "--private_mutations",
                        help="Add this flag to search for private mutations. These are mutations that are not part of"
                             " a haplogroup but are annotated in dbsnp or differentiate from the reference genome.",
                        action="store_true")

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


def main_fastq(
    args: argparse.Namespace,
    out_folder: Path
):
    files = get_files_with_extension(args.fastq, '.fastq')
    reference = get_reference_path(args.reference_genome, True)
    bam_folder = out_folder / "bam_files"
    try:
        os.mkdir(bam_folder)
    except IOError:
        pass
    LOG.info("Creating bam files from fastq files...")
    for fastq_file in files:
        LOG.info(f"Starting with running for {fastq_file}")
        sam_file = bam_folder / "temp_fastq_sam.sam"

        fastq_cmd = f"minimap2 -ax sr -t {args.threads} {reference} {fastq_file} > {sam_file}"
        call_command(fastq_cmd)
        bam_file = bam_folder / (fastq_file.name.rsplit(".", 1)[0] + ".bam")
        cmd = "samtools view -@ {} -bS {} | samtools sort -@ {} -m 2G -o {}".format(args.threads, sam_file,
                                                                                    args.threads, bam_file)
        call_command(cmd)
        cmd = "samtools index -@ {} {}".format(args.threads, bam_file)
        call_command(cmd)
        os.remove(sam_file)
    args.bamfile = bam_folder
    main_bam_cram(args, out_folder, True)


def main_bam_cram(
    args: argparse.Namespace,
    base_out_folder: Path,
    is_bam: bool
):
    files = get_files_with_extension(args.bamfile, '.bam')
    # make sure the reference genome is present before doing something else
    check_reference(args.reference_genome)
    with multiprocessing.Pool(processes=args.threads) as p:
        p.map(partial(run_bam_cram, args, base_out_folder, is_bam), files)


def run_bam_cram(
    args: argparse.Namespace,
    base_out_folder: Path,
    is_bam: bool,
    input_file: Path
):
    LOG.info(f"Starting with running for {input_file}")
    output_dir = base_out_folder / input_file.name.split(".")[0]
    safe_create_dir(output_dir, args.force)
    general_info_list = samtools(output_dir, input_file, is_bam, args)
    write_info_file(output_dir, general_info_list)
    if args.private_mutations:
        find_private_mutations(output_dir, input_file, args, is_bam)
    LOG.info(f"Finished running for {input_file.name}")
    print()


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


def check_reference(
    requested_version: str,
):
    reference_file = get_reference_path(requested_version, True)
    if not reference_file.exists():
        LOG.info(f"No reference genome version was found. Downloading the {requested_version} reference genome. This "
                 f"should be a one time thing.")
        download_reference.main(requested_version)
        LOG.info("Finished downloading the reference genome.")


def get_reference_path(
    requested_version: str,
    is_full: bool
) -> Union[Path, None]:
    if is_full:
        reference_file = yleaf_constants.DATA_FOLDER / requested_version / yleaf_constants.FULL_REF_FILE
    else:
        reference_file = yleaf_constants.DATA_FOLDER / requested_version / yleaf_constants.Y_REF_FILE

    return reference_file


def samtools(
    output_folder: Path,
    path_file: Path,
    is_bam_pathfile: bool,
    args: argparse.Namespace,
) -> List[str]:

    outputfile = output_folder / (output_folder.name + ".out")
    fmf_output = output_folder / (output_folder.name + ".fmf")
    pileupfile = output_folder / "temp_haplogroup_pileup.pu"
    reference = get_reference_path(args.reference_genome, True) if not is_bam_pathfile else None

    if is_bam_pathfile:
        if not any([Path(str(path_file) + ".bai").exists(), Path(str(path_file).rstrip(".bam") + '.bai').exists()]):
            cmd = "samtools index -@{} {}".format(args.threads, path_file)
            call_command(cmd)
    else:
        if not any([Path(str(path_file) + ".crai").exists(), Path(str(path_file).rstrip(".cram") + '.crai').exists()]):
            cmd = "samtools index -@{} {}".format(args.threads, path_file)
            call_command(cmd)
    header, mapped, unmapped = chromosome_table(path_file, output_folder, output_folder.name)

    position_file = get_position_file(args.reference_genome, args.use_old)

    bed = output_folder / "temp_position_bed.bed"
    write_bed_file(bed, position_file, header)

    execute_mpileup(bed, path_file, pileupfile, args.quality_thresh, reference)

    general_info_list = ["Total of mapped reads: " + str(mapped), "Total of unmapped reads: " + str(unmapped)]

    extract_haplogroups(position_file, args.reads_treshold, args.base_majority,
                        pileupfile, fmf_output, outputfile, is_bam_pathfile, args.use_old, general_info_list)

    os.remove(pileupfile)
    os.remove(bed)

    LOG.info("Finished extracting haplogroups")
    return general_info_list


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


def get_position_file(
    reference_name: str,
    use_old: bool
) -> Path:
    if use_old:
        position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.NEW_POSITION_FILE
    else:
        position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.NEW_POSITION_FILE
    return position_file


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
    cmd += f" -AQ{quality_thresh}q1 {str(bam_file)} > {str(pileupfile)}"
    call_command(cmd)


def extract_haplogroups(
    path_markerfile: Path,
    reads_thresh: float,
    base_majority: int,
    path_pileupfile: Path,
    fmf_output: Path,
    outputfile: Path,
    is_bam_file: bool,
    use_old: bool,
    general_info_list: List[str]
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
    general_info_list.append("Valid markers: " + str(markerfile_len))

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
    called_base = df_freq_table.columns[list_col_indices]  # noqa
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

    general_info_list.append("Markers with zero reads: " + str(len(df_belowzero)))
    general_info_list.append(
        "Markers below the read threshold {" + str(reads_thresh) + "}: " + str(len(df_readsthreshold)))
    general_info_list.append(
        "Markers below the base majority threshold {" + str(base_majority) + "}: " + str(len(df_basemajority)))
    general_info_list.append("Markers with discordant genotype: " + str(len(df_discordantgenotype)))
    general_info_list.append("Markers without haplogroup information: " + str(len(df_fmf)))
    general_info_list.append("Markers with haplogroup information: " + str(len(df_out)))

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

    tree = Tree(yleaf_constants.DATA_FOLDER / yleaf_constants.HG_PREDICTION_FOLDER / yleaf_constants.TREE_FILE)
    with open(outputfile, "w") as f:
        f.write('\t'.join(["chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads",
                           "called_perc", "called_base", "state", "depth\n"]))
        for node_key in tree.node_mapping:
            if node_key not in mappable_df:
                continue
            depth = tree.get(node_key).depth
            for lst in mappable_df[node_key]:
                f.write('\t'.join(map(str, lst)) + f"\t{depth}\n")


def replace_with_bases(
    base: str,
    read_result: str
) -> str:
    return read_result.replace(",", base[0]).replace(".", base[0])


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
    fastadict = {"A": 0, "T": 0, "G": 0, "C": 0, "-": 0, "+": 0, "*": 0}
    sequence = sequence.upper()
    index = 0
    while index < len(sequence):
        char = sequence[index]
        if char in fastadict:
            fastadict[char] += 1
            index += 1
        elif char == "^":
            index += 2
        elif char in {"-", "+"}:
            index += 1
            digit, index = find_digit(sequence, index)
            index += digit
        else:
            index += 1
    fastadict["-"] += fastadict["*"]
    del fastadict["*"]
    return fastadict


def find_digit(
    sequence: str,
    index: int
) -> Tuple[int, int]:
    # first is always a digit
    nr = [sequence[index]]
    index += 1
    while True:
        char = sequence[index]
        # this seems to be faster than isdigit()
        if char in NUM_SET:
            nr.append(char)
            index += 1
            continue
        return int(''.join(nr)), index


def write_info_file(
    folder: Path,
    general_info_list: List[str]
):
    try:
        with open(folder / (folder.name + ".info"), "a") as f:
            for marker in general_info_list:
                f.write(marker)
                f.write("\n")
    except IOError:
        LOG.warning("Failed to write .info file")


def find_private_mutations(
    output_folder: Path,
    path_file: Path,
    args: argparse.Namespace,
    is_bam: bool
):
    # identify mutations not part of haplogroups that are annotated in dbsnp or differ from the reference genome
    LOG.info("Starting with extracting private mutations...")
    full_reference = get_reference_path(args.reference_genome, True) if not is_bam else None
    snp_reference_file = get_reference_path(args.reference_genome, False)
    snp_database_file = yleaf_constants.DATA_FOLDER / args.reference_genome / yleaf_constants.SNP_DATA_FILE

    # run mpileup
    pileup_file = output_folder / "temp_private_mutation_pileup.pu"
    execute_mpileup(None, path_file, pileup_file, args.quality_thresh, full_reference)

    LOG.info("Loading reference files")

    position_file = get_position_file(args.reference_genome, args.use_old)
    filter_positions = load_filter_data(position_file)
    known_snps = load_snp_database_file(snp_database_file, args.maximum_mutation_rate)
    ychrom_reference = load_reference_file(snp_reference_file)

    LOG.info("Finding private mutations...")
    private_mutations = []
    confirmed_private_mutations = []
    with open(pileup_file) as f:
        for index, line in enumerate(f):
            try:
                chrom, position, ref_base, count, aligned, quality = line.strip().split()
            except ValueError:
                LOG.warning(f"failed to read line {index} of pileupfile")
                continue
            if chrom != "chrY":
                continue
            # not enough reads
            if int(count) < args.reads_treshold:
                continue
            if not is_bam:
                aligned = replace_with_bases(ref_base, aligned)
            if position in filter_positions:
                continue

            # not annotated in dbsnp
            if position not in known_snps and snp_reference_file is not None:
                freq_dict = get_frequencies(aligned)
                actual_allele, allele_count = max(freq_dict.items(), key=lambda x: x[1])
                called_percentage = round(allele_count / sum(freq_dict.values()) * 100, 2)
                # not a high enough base majority measured, cannot be sure of real allele
                if called_percentage < args.base_majority:
                    continue
                ref_base = ychrom_reference[int(position) - 1]

                # do not match against insertions or repeat regions (lower case)
                if ref_base not in ACCEPTED_REF_BASES or actual_allele not in ACCEPTED_REF_BASES:
                    continue
                if ref_base == actual_allele:
                    continue
                private_mutations.append(f"{chrom}\t{position}\t-\t{ref_base}->{actual_allele}\t{ref_base}\t"
                                         f"{actual_allele}\t{allele_count}\t{called_percentage}\tNA\n")
            elif position in known_snps and snp_database_file is not None:
                freq_dict = get_frequencies(aligned)
                actual_allele, allele_count = max(freq_dict.items(), key=lambda x: x[1])
                called_percentage = round(allele_count / sum(freq_dict.values()) * 100, 2)
                # not a high enough base majority measured, cannot be sure of real allele
                if called_percentage < args.base_majority:
                    continue
                possible_minor_alleles = {dct["minor_allele"] for dct in known_snps[position]}
                # measured allele is not the dbsnp allele
                if actual_allele not in possible_minor_alleles:
                    continue

                matched_pos_dct = [dct for dct in known_snps[position]
                                   if dct["minor_allele"] == actual_allele][0]
                rs, major_allele, minor_allele, frequency = matched_pos_dct.values()
                confirmed_private_mutations.append(f"{chrom}\t{position}\t{rs}\t{major_allele}->{actual_allele}\t"
                                                   f"{major_allele}\t{actual_allele}\t{allele_count}\t"
                                                   f"{called_percentage}\t{frequency}\n")
    os.remove(pileup_file)
    with open(output_folder / "private_mutations.csv", "w") as f:
        f.write(f"chrom\tposition\trs_nr\tmutation\treference\tactual\treads\tcalled_percentage\tfrequency\n")
        f.write(''.join(confirmed_private_mutations))
        f.write(''.join(private_mutations))
    LOG.info("Finished extracting private mutations")


def load_filter_data(
    path: Path
) -> Set[str]:
    global CACHED_POSITION_DATA
    if CACHED_POSITION_DATA is None:
        CACHED_POSITION_DATA = set()
        with open(path) as f:
            for line in f:
                CACHED_POSITION_DATA.add(line.strip().split("\t")[3])
    return CACHED_POSITION_DATA


def load_snp_database_file(
    path: Path,
    maximum_mutation_rate: float
) -> Union[Dict[str, List[Dict[str, str]]], None]:
    global CACHED_SNP_DATABASE

    if path is not None:
        if CACHED_SNP_DATABASE is None:
            CACHED_SNP_DATABASE = defaultdict(list)
            with open(path) as f:
                f.readline()
                for line in f:
                    rs, position, major_allele, minor_allele, frequency = line.strip().split(",")
                    if float(frequency) > maximum_mutation_rate:
                        continue
                    CACHED_SNP_DATABASE[position].append({"rs": rs, "major_allele": major_allele,
                                                          "minor_allele": minor_allele, "frequency": frequency})
            CACHED_SNP_DATABASE = dict(CACHED_SNP_DATABASE)
        return CACHED_SNP_DATABASE
    else:
        return None


def load_reference_file(
    path: Path
) -> Union[List[str], None]:
    global CACHED_REFERENCE_FILE
    if path is not None:
        if CACHED_REFERENCE_FILE is None:
            CACHED_REFERENCE_FILE = []
            with open(path) as f:
                for line in f:
                    if line.startswith(">"):
                        continue
                    CACHED_REFERENCE_FILE.extend(line.strip())
        return CACHED_REFERENCE_FILE
    else:
        return None


def predict_haplogroup(
    path_file: Path,
    output: Path,
    use_old: bool
):
    if use_old:
        script = yleaf_constants.SRC_FOLDER / "old_predict_haplogroup.py"
        cmd = "python {} -i {} -o {}".format(script, path_file, output)
        call_command(cmd)
    else:
        from yleaf import predict_haplogroup
        namespace = argparse.Namespace(input=path_file, outfile=output,
                                       minimum_score=predict_haplogroup.DEFAULT_MIN_SCORE)
        predict_haplogroup.main(namespace)


def draw_haplogroups(
    haplogroup_file: Path,
    collapsed_draw_mode: bool
):
    # make sure that it is only imported if requested by user
    from yleaf import haplogroup_tree_image
    namespace = argparse.Namespace(input=haplogroup_file, collapse_mode=collapsed_draw_mode,
                                   outfile=haplogroup_file.parent / HAPLOGROUP_IMAGE_FILE_NAME)
    haplogroup_tree_image.main(namespace)


if __name__ == "__main__":
    main()