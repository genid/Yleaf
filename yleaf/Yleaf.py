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
import time
import datetime

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


class MyFormatter(logging.Formatter):
    """
    Copied from MultiGeneBlast (my code)
    """

    def __init__(self, fmt, starttime=time.time()):
        logging.Formatter.__init__(self, fmt)
        self._start_time = starttime

    def format(self, record):
        """
        Overwrite of the format function that prints the passed time and adds
        current time to the existing format
        :See: logging.Formatter.format()
        """
        record.passedTime = "{:.3f}".format(time.time() - self._start_time)
        record.currentTime = datetime.datetime.now().time()
        return super(MyFormatter, self).format(record)


def run_vcf(
        path_markerfile: Path,
        base_out_folder: Path,
        args: argparse.Namespace,
        sample_vcf_file: Path,
        out_suffix: str = '',
        tree_name: str = None,
):
    LOG.debug("Starting with extracting haplogroups...")
    markerfile = pd.read_csv(path_markerfile, header=None, sep="\t",
                             names=["chr", "marker_name", "haplogroup", "pos", "mutation", "anc",
                                    "der"],
                             dtype={"chr": str,
                                    "marker_name": str,
                                    "haplogroup": str,
                                    "pos": int,
                                    "mutation": str,
                                    "anc": str,
                                    "der": str,
                                    }).drop_duplicates(subset='pos', keep='first', inplace=False)

    sample_vcf_folder = base_out_folder / (sample_vcf_file.name.replace(".vcf.gz", ""))
    safe_create_dir(sample_vcf_folder, args.force if not out_suffix else False, reuse_pileup=bool(out_suffix))

    sample_vcf_file_txt = sample_vcf_folder / (sample_vcf_file.name.replace(".vcf.gz", ".txt"))
    cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' {sample_vcf_file} > {sample_vcf_file_txt}"
    call_command(cmd)

    pileupfile = pd.read_csv(sample_vcf_file_txt,
                             dtype=str, header=None, sep="\t")

    # remove sample_vcf_file_txt
    cmd = f"rm {sample_vcf_file_txt}"
    call_command(cmd)

    pileupfile.columns = ['chr', 'pos', 'refbase', 'altbase', 'reads']
    pileupfile['pos'] = pileupfile['pos'].astype(int)

    pileupfile['altbase'] = pileupfile['altbase'].str.split(',')
    pileupfile['reads'] = pileupfile['reads'].str.split(',')

    pileupfile['ref_reads'] = pileupfile['reads'].apply(lambda x: x[0])
    pileupfile['alt_reads'] = pileupfile['reads'].apply(lambda x: x[1:])

    pileupfile['alt_reads_dict'] = pileupfile.apply(lambda row: dict(zip(row['altbase'], row['alt_reads'])), axis=1)
    pileupfile['alt_reads_dict'] = pileupfile['alt_reads_dict'].apply(lambda x: {k: int(v) for k, v in x.items()})
    pileupfile['highest_alt_reads'] = pileupfile['alt_reads_dict'].apply(lambda x: max(x.values()) if len(x) > 0 else 0)
    pileupfile['highest_alt_reads_base'] = pileupfile['alt_reads_dict'].apply(
        lambda x: max(x, key=x.get) if len(x) > 0 else 'NA')
    pileupfile['total_reads'] = pileupfile.apply(lambda row: int(row['ref_reads']) + row['highest_alt_reads'], axis=1)
    pileupfile['called_ref_perc'] = pileupfile.apply(
        lambda row: round((int(row['ref_reads']) / row['total_reads']) * 100, 1) if row['total_reads'] > 0 else 0,
        axis=1)
    pileupfile['called_alt_perc'] = pileupfile.apply(
        lambda row: round((row['highest_alt_reads'] / row['total_reads']) * 100, 1) if row['total_reads'] > 0 else 0,
        axis=1)

    pileupfile['called_base'] = pileupfile.apply(
        lambda row: row['refbase'] if row['called_ref_perc'] >= row['called_alt_perc'] else row[
            'highest_alt_reads_base'], axis=1)
    pileupfile['called_perc'] = pileupfile.apply(
        lambda row: row['called_ref_perc'] if row['called_ref_perc'] >= row['called_alt_perc'] else row[
            'called_alt_perc'], axis=1).astype(float)
    pileupfile['called_reads'] = pileupfile.apply(
        lambda row: row['ref_reads'] if row['called_ref_perc'] >= row['called_alt_perc'] else row['highest_alt_reads'],
        axis=1).astype(int)

    # Save full marker list before intersection — needed to infer states for positions
    # absent from the VCF (where the sample matches the reference genome).
    full_markerfile = markerfile

    intersect_pos = np.intersect1d(pileupfile['pos'], markerfile['pos'])
    markerfile = markerfile.loc[markerfile['pos'].isin(intersect_pos)]
    markerfile = markerfile.sort_values(by=['pos'])
    pileupfile = pileupfile.loc[pileupfile['pos'].isin(intersect_pos)]

    pileupfile = pileupfile.drop(['chr'], axis=1)
    df = pd.merge(markerfile, pileupfile, on='pos')

    df['state'] = df.apply(
        lambda row: 'A' if row['called_base'] == row['anc'] else 'D' if row['called_base'] == row['der'] else 'NA',
        axis=1)
    df['bool_state'] = df.apply(
        lambda row: True if (row['called_base'] == row['anc'] or row['called_base'] == row['der']) else False, axis=1)

    markerfile_len = len(markerfile)

    # valid markers from positionsfile.txt
    general_info_list = ["Total of mapped reads: VCF", "Total of unmapped reads: VCF"]
    general_info_list += ["Valid markers: " + str(markerfile_len)]

    index_belowzero = df[df["called_reads"] == 0].index
    df_belowzero = df[df.index.isin(index_belowzero)]
    df_belowzero = df_belowzero.drop(['refbase', 'altbase'], axis=1)
    df_belowzero["called_perc"] = "NA"
    df_belowzero["called_base"] = "NA"
    df_belowzero["state"] = "NA"
    df_belowzero["Description"] = "Position with zero reads"

    df = df[~df.index.isin(index_belowzero)]
    bool_list_state = df[df["bool_state"] == False].index

    # discordant genotypes
    df_discordantgenotype = df[df.index.isin(bool_list_state)]
    df_discordantgenotype = df_discordantgenotype.drop(["bool_state"], axis=1)
    df_discordantgenotype["state"] = "NA"
    df_discordantgenotype["Description"] = "Discordant genotype"
    df = df[~df.index.isin(bool_list_state)]

    reads_thresh = int(args.reads_treshold)

    # read threshold
    df_readsthreshold = df[df["called_reads"] < reads_thresh]
    df_readsthreshold["Description"] = "Below read threshold"
    df = df[df["called_reads"] >= reads_thresh]

    base_majority = float(args.base_majority)

    # filter by base percentage
    df_basemajority = df[df["called_perc"] < base_majority]
    df_basemajority["Description"] = "Below base majority"
    df = df[df["called_perc"] >= base_majority]

    df_fmf = pd.concat([df_belowzero, df_readsthreshold, df_basemajority, df_discordantgenotype], axis=0, sort=True)
    df_fmf['reads'] = df_fmf['called_reads']
    df_fmf = df_fmf[['chr', 'pos', 'marker_name', 'haplogroup', 'mutation', 'anc', 'der', 'reads',
                     'called_perc', 'called_base', 'state', 'Description']]

    df_out = df
    df_out['reads'] = df_out['called_reads']
    df_out = df_out[['chr', 'pos', 'marker_name', 'haplogroup', 'mutation', 'anc', 'der', 'reads',
                     'called_perc', 'called_base', 'state']]

    # Infer states for marker positions absent from the VCF.
    # Standard VCFs only contain variant positions; ancestral positions where the sample
    # matches the reference are not emitted. By looking up the reference allele at those
    # positions from the downloaded chrY FASTA we can recover the vast majority of markers.
    ref_inferred_count = 0
    chry_ref_path = get_reference_path(args.reference_genome, False)
    if chry_ref_path and chry_ref_path.exists() and os.path.getsize(chry_ref_path) > 0:
        try:
            with open(chry_ref_path) as _ref_f:
                chry_seq = ''.join(line.strip().upper() for line in _ref_f if not line.startswith('>'))
            covered_pos_set = set(intersect_pos)
            missing_markers = full_markerfile[~full_markerfile['pos'].isin(covered_pos_set)].copy()
            if len(missing_markers) > 0:
                missing_markers['ref_allele'] = missing_markers['pos'].apply(
                    lambda p: chry_seq[p - 1] if 0 < p <= len(chry_seq) else 'N')
                missing_markers = missing_markers[missing_markers.apply(
                    lambda row: row['ref_allele'] in ACCEPTED_REF_BASES
                                and row['ref_allele'] in (row['anc'], row['der']), axis=1)].copy()
                if len(missing_markers) > 0:
                    missing_markers['state'] = missing_markers.apply(
                        lambda row: 'A' if row['ref_allele'] == row['anc'] else 'D', axis=1)
                    missing_markers['reads'] = 0
                    missing_markers['called_perc'] = 100.0
                    missing_markers['called_base'] = missing_markers['ref_allele']
                    ref_inferred_df = missing_markers[[
                        'chr', 'pos', 'marker_name', 'haplogroup', 'mutation', 'anc', 'der',
                        'reads', 'called_perc', 'called_base', 'state']]
                    df_out = pd.concat([df_out, ref_inferred_df], axis=0, sort=False)
                    ref_inferred_count = len(ref_inferred_df)
                    LOG.debug(
                        f"Reference inference: added {ref_inferred_count} markers "
                        f"({(missing_markers['state'] == 'A').sum()} A, "
                        f"{(missing_markers['state'] == 'D').sum()} D)")
        except Exception as e:
            LOG.warning(f"Could not infer states from reference genome for VCF input: {e}")
    else:
        LOG.warning(
            "chrY reference FASTA not available — VCF prediction will use only variant positions. "
            "Run Yleaf once with -bam to download the reference, or use -bam/-cram input.")

    general_info_list.append("Markers with zero reads: " + str(len(df_belowzero)))
    general_info_list.append(
        "Markers below the read threshold {" + str(reads_thresh) + "}: " + str(len(df_readsthreshold)))
    general_info_list.append(
        "Markers below the base majority threshold {" + str(base_majority) + "}: " + str(len(df_basemajority)))
    general_info_list.append("Markers with discordant genotype: " + str(len(df_discordantgenotype)))
    general_info_list.append("Markers without haplogroup information: " + str(len(df_fmf)))
    general_info_list.append("Reference-inferred markers: " + str(ref_inferred_count))
    general_info_list.append("Markers with haplogroup information: " + str(len(df_out)))

    write_info_file(sample_vcf_folder, general_info_list, suffix=f"{out_suffix}.info" if out_suffix else ".info")

    use_old = args.use_old
    outputfile = sample_vcf_folder / (sample_vcf_file.name.replace(".vcf.gz", f"{out_suffix}.out"))
    fmf_output = sample_vcf_folder / (sample_vcf_file.name.replace(".vcf.gz", f"{out_suffix}.fmf"))

    if use_old:
        df_out = df_out.sort_values(by=['haplogroup'], ascending=True)
        df_out = df_out[
            ["chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads", "called_perc",
             "called_base",
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

    tree = Tree(get_tree_path(tree_name if tree_name else args.tree))
    with open(outputfile, "w") as f:
        f.write('\t'.join(["chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads",
                           "called_perc", "called_base", "state", "depth\n"]))
        for node_key in tree.node_mapping:
            if node_key not in mappable_df:
                continue
            depth = tree.get(node_key).depth
            for lst in mappable_df[node_key]:
                f.write('\t'.join(map(str, lst)) + f"\t{depth}\n")

    LOG.info(f"Finished extracting genotypes for {sample_vcf_file.name}")


# ── PLINK / SNP array support ─────────────────────────────────────────────────

_PLINK_Y_CHROM_CODES: Set[str] = {'24', 'Y', 'chrY', 'y', 'chry', 'ChrY'}


def _plink_detect_base(plink_input: Path) -> Path:
    """Return the base path (without extension) for a PLINK file set."""
    if plink_input.suffix.lower() in ('.ped', '.bed', '.bim', '.fam', '.map'):
        return plink_input.with_suffix('')
    return plink_input  # already a base name


def _parse_plink_bim(bim_path: Path):
    """Parse .bim file. Returns (y_snps, total_snp_count).
    y_snps: list of (snp_index, bp_pos, allele1, allele2) for Y-chr entries."""
    y_snps = []
    snp_idx = 0
    with open(bim_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6 and parts[0] in _PLINK_Y_CHROM_CODES:
                try:
                    y_snps.append((snp_idx, int(parts[3]), parts[4], parts[5]))
                except ValueError:
                    pass
            snp_idx += 1
    return y_snps, snp_idx


def _parse_plink_fam(fam_path: Path) -> List[str]:
    """Parse .fam file. Returns list of sample IIDs."""
    samples = []
    with open(fam_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                samples.append(parts[1])
    return samples


def _parse_plink_binary(base: Path) -> Dict[str, Dict[int, Union[str, None]]]:
    """Parse .bed/.bim/.fam. Returns {sample_id: {pos: called_base}} for Y-chr SNPs.
    Heterozygous or missing calls are omitted (None not stored)."""
    y_snps, total_snps = _parse_plink_bim(base.with_suffix('.bim'))
    samples = _parse_plink_fam(base.with_suffix('.fam'))
    result: Dict[str, Dict[int, Union[str, None]]] = {s: {} for s in samples}
    if not y_snps or not samples:
        return result

    n_samples = len(samples)
    bytes_per_snp = (n_samples + 3) // 4

    with open(base.with_suffix('.bed'), 'rb') as f:
        magic = f.read(3)
        if len(magic) < 3 or magic[:2] != b'\x6c\x1b':
            raise ValueError(f"Not a valid PLINK .bed file: {base}.bed")
        if magic[2] != 0x01:
            raise ValueError("Only SNP-major PLINK .bed format is supported (byte 3 must be 0x01)")
        bed_data = f.read()

    for snp_idx, pos, a1, a2 in y_snps:
        byte_start = snp_idx * bytes_per_snp
        snp_bytes = bed_data[byte_start:byte_start + bytes_per_snp]
        for sample_idx, sample_id in enumerate(samples):
            byte_val = snp_bytes[sample_idx // 4]
            geno_bits = (byte_val >> ((sample_idx % 4) * 2)) & 0b11
            # 00 = hom A1 / 01 = missing / 10 = het / 11 = hom A2
            if geno_bits == 0b00 and a1 not in ('0', '.'):
                result[sample_id][pos] = a1
            elif geno_bits == 0b11 and a2 not in ('0', '.'):
                result[sample_id][pos] = a2
            # missing (01) and het (10) → not stored

    return result


def _parse_plink_text(base: Path) -> Dict[str, Dict[int, Union[str, None]]]:
    """Parse .ped/.map. Returns {sample_id: {pos: called_base}} for Y-chr SNPs.
    Heterozygous or missing calls are omitted."""
    y_snp_indices: List[Tuple[int, int]] = []  # (col_index, bp_pos)
    snp_idx = 0
    with open(base.with_suffix('.map')) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4 and parts[0] in _PLINK_Y_CHROM_CODES:
                try:
                    y_snp_indices.append((snp_idx, int(parts[3])))
                except ValueError:
                    pass
            snp_idx += 1

    result: Dict[str, Dict[int, Union[str, None]]] = {}
    if not y_snp_indices:
        return result

    with open(base.with_suffix('.ped')) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 6:
                continue
            sample_id = parts[1]
            genos = parts[6:]
            result[sample_id] = {}
            for snp_idx, pos in y_snp_indices:
                col = snp_idx * 2
                if col + 1 >= len(genos):
                    continue
                a1, a2 = genos[col], genos[col + 1]
                if a1 == '0' or a2 == '0':
                    continue  # missing
                if a1 == a2:
                    result[sample_id][pos] = a1  # haploid Y call
                # heterozygous on Y → ambiguous, skip

    return result


def _write_plink_sample_out(
        sample_id: str,
        y_genos: Dict[int, str],
        markerfile: pd.DataFrame,
        sample_folder: Path,
        args: argparse.Namespace,
        tree_name: str,
        out_suffix: str = '',
):
    """Write .out and .info files for one sample from PLINK genotypes."""
    tree = Tree(get_tree_path(tree_name))

    rows_by_hg: Dict[str, list] = {}
    n_called = n_missing = n_discordant = 0

    for row in markerfile.itertuples(index=False):
        called = y_genos.get(row.pos)
        if called is None:
            n_missing += 1
            continue
        if called == row.anc:
            state = 'A'
        elif called == row.der:
            state = 'D'
        else:
            n_discordant += 1
            continue
        hg = row.haplogroup
        if hg not in rows_by_hg:
            rows_by_hg[hg] = []
        rows_by_hg[hg].append([row.chr, row.pos, row.marker_name, hg,
                                row.mutation, row.anc, row.der,
                                1, 100.0, called, state])
        n_called += 1

    out_path = sample_folder / f"{sample_id}{out_suffix}.out"
    header_cols = ["chr", "pos", "marker_name", "haplogroup", "mutation",
                   "anc", "der", "reads", "called_perc", "called_base", "state"]

    if args.use_old:
        with open(out_path, 'w') as f:
            f.write('\t'.join(header_cols) + '\n')
            for hg in sorted(rows_by_hg):
                for r in rows_by_hg[hg]:
                    f.write('\t'.join(map(str, r)) + '\n')
    else:
        with open(out_path, 'w') as f:
            f.write('\t'.join(header_cols + ["depth"]) + '\n')
            for node_key in tree.node_mapping:
                if node_key not in rows_by_hg:
                    continue
                depth = tree.get(node_key).depth
                for r in rows_by_hg[node_key]:
                    f.write('\t'.join(map(str, r)) + f'\t{depth}\n')

    info_suffix = f"{out_suffix}.info" if out_suffix else ".info"
    write_info_file(sample_folder, [
        "Total of mapped reads: PLINK array",
        "Total of unmapped reads: PLINK array",
        f"Valid markers: {len(markerfile)}",
        f"Markers not on array: {n_missing}",
        f"Markers with discordant genotype: {n_discordant}",
        f"Markers without haplogroup information: {n_missing + n_discordant}",
        f"Markers with haplogroup information: {n_called}",
    ], suffix=info_suffix)

    LOG.info(f"  {sample_id}: {n_called} markers called, "
             f"{n_missing} not on array, {n_discordant} discordant")


def _load_plink_genotypes(args: argparse.Namespace) -> Dict[str, Dict[int, str]]:
    """Parse PLINK file and return Y genotypes. Shared by single- and multi-tree paths."""
    base = _plink_detect_base(Path(args.plinkfile))
    is_binary = base.with_suffix('.bed').exists()
    fmt = 'binary (.bed/.bim/.fam)' if is_binary else 'text (.ped/.map)'
    LOG.info(f"Parsing PLINK {fmt}: {base}")
    y_genos = _parse_plink_binary(base) if is_binary else _parse_plink_text(base)
    if not y_genos:
        LOG.warning("No Y-chromosome SNPs found. Check chromosome encoding "
                    "(expected: 24, Y, or chrY in column 1 of .bim/.map).")
    else:
        LOG.info(f"Found {len(y_genos)} sample(s) in PLINK file")
    return y_genos


def main_plink(args: argparse.Namespace, base_out_folder: Path):
    """Single-tree PLINK entry point."""
    y_genos = _load_plink_genotypes(args)
    if not y_genos:
        return
    markerfile = pd.read_csv(
        get_position_file(args.reference_genome, args.use_old, args.ancient_DNA, args.tree),
        header=None, sep="\t",
        names=["chr", "marker_name", "haplogroup", "pos", "mutation", "anc", "der"],
        dtype={"chr": str, "marker_name": str, "haplogroup": str,
               "pos": int, "mutation": str, "anc": str, "der": str},
    ).drop_duplicates(subset='pos', keep='first')

    for sample_id, sample_y_genos in y_genos.items():
        sample_folder = base_out_folder / sample_id
        safe_create_dir(sample_folder, args.force)
        _write_plink_sample_out(sample_id, sample_y_genos, markerfile,
                                sample_folder, args, tree_name=args.tree)


def main_plink_multi_tree(args: argparse.Namespace, base_out_folder: Path, trees: List[str]):
    """Multi-tree PLINK entry point. Parses the PLINK file once, predicts per tree."""
    y_genos = _load_plink_genotypes(args)
    if not y_genos:
        return

    # Create all sample folders upfront (avoids force-deletion on second tree pass)
    for sample_id in y_genos:
        safe_create_dir(base_out_folder / sample_id, args.force)

    for tree in trees:
        LOG.info(f"Extracting markers for tree: {tree}")
        markerfile = pd.read_csv(
            get_position_file(args.reference_genome, args.use_old, args.ancient_DNA, tree),
            header=None, sep="\t",
            names=["chr", "marker_name", "haplogroup", "pos", "mutation", "anc", "der"],
            dtype={"chr": str, "marker_name": str, "haplogroup": str,
                   "pos": int, "mutation": str, "anc": str, "der": str},
        ).drop_duplicates(subset='pos', keep='first')
        for sample_id, sample_y_genos in y_genos.items():
            _write_plink_sample_out(sample_id, sample_y_genos, markerfile,
                                    base_out_folder / sample_id, args,
                                    tree_name=tree, out_suffix=f'.{tree}')

    combined_out = base_out_folder / "hg_prediction_combined.hg"
    write_combined_prediction_table(
        base_out_folder, trees, combined_out,
        args.use_old, args.prediction_quality, args.threads,
        draw_hg=args.draw_haplogroups,
        collapsed_draw_mode=args.collapsed_draw_mode,
    )


def main_vcf_split(
        position_bed_file: Path,
        base_out_folder: Path,
        args: argparse.Namespace,
        vcf_file: Path
):
    # first sort the vcf file
    sorted_vcf_file = base_out_folder / (vcf_file.name.replace(".vcf.gz", ".sorted.vcf.gz"))
    cmd = f"bcftools sort -O z -o {sorted_vcf_file} {vcf_file}"
    try:
        call_command(cmd)
    except SystemExit:
        LOG.error(f"Failed to sort the vcf file {vcf_file.name}. Skipping...")
        return None

    # next index the vcf file
    cmd = f"bcftools index -f {sorted_vcf_file}"
    call_command(cmd)

    # get chromosome annotation using bcftools query -f '%CHROM\n' sorted.vcf.gz | uniq
    cmd = f"bcftools query -f '%CHROM\n' {sorted_vcf_file} | uniq"
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        LOG.error(f"Call: '{cmd}' failed. Reason given: '{stderr.decode('utf-8')}'")
        raise SystemExit("Failed command execution")
    chromosomes = stdout.decode('utf-8').strip().split("\n")
    chry = [x for x in chromosomes if "y" in x.lower()]
    if len(chry) == 0:
        LOG.error("Unable to find Y-chromosome in the vcf file. Exiting...")
        raise SystemExit("Unable to find Y-chromosome in the vcf file.")
    elif len(chry) > 1:
        LOG.error("Multiple Y-chromosome annotations found in the vcf file. Exiting...")
        raise SystemExit("Multiple Y-chromosome annotations found in the vcf file.")
    else:
        # make new position_bed_file with correct chrY annotation
        new_position_bed_file = base_out_folder / (vcf_file.name.replace(".vcf.gz", "temp_position_bed.bed"))
        with open(position_bed_file, "r") as f:
            with open(new_position_bed_file, "w") as f2:
                for line in f:
                    line = line.replace("chrY", chry[0])
                    f2.write(line)

    # filter the vcf file using the reference bed file
    filtered_vcf_file = base_out_folder / "filtered_vcf_files" / (
        sorted_vcf_file.name.replace(".sorted.vcf.gz", ".filtered.vcf.gz"))
    cmd = f"bcftools view -O z -R {new_position_bed_file} {sorted_vcf_file} > {filtered_vcf_file}"
    call_command(cmd)

    # remover temp_position_bed.bed
    cmd = f"rm {new_position_bed_file}"
    call_command(cmd)

    # remove sorted.vcf.gz and sorted.vcf.gz.csi
    cmd = f"rm {sorted_vcf_file}"
    call_command(cmd)
    cmd = f"rm {sorted_vcf_file}.csi"
    call_command(cmd)

    # check number of samples in the vcf file
    cmd = f"bcftools query -l {filtered_vcf_file} | wc -l"
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        LOG.error(f"Call: '{cmd}' failed. Reason given: '{stderr.decode('utf-8')}'")
        raise SystemExit("Failed command execution")
    num_samples = int(stdout.decode('utf-8').strip())

    if num_samples > 1:
        # split the vcf file into separate files for each sample
        split_vcf_folder = base_out_folder / (vcf_file.name.replace(".vcf.gz", "_split"))
        safe_create_dir(split_vcf_folder, args.force)
        cmd = f"bcftools +split {filtered_vcf_file} -Oz -o {split_vcf_folder}"
        call_command(cmd)
        sample_vcf_files = get_files_with_extension(split_vcf_folder, '.vcf.gz')
    elif num_samples == 1:
        sample_vcf_files = [filtered_vcf_file]
    else:
        LOG.error("No samples found in the vcf file. Exiting...")
        raise SystemExit("No samples found in the vcf file.")

    return sample_vcf_files


def _write_merged_vcf_bed(merged_bed: Path, trees: List[str], reference_genome: str,
                          use_old: bool, ancient_dna: bool) -> None:
    """Merge BED files from multiple trees into a union BED for VCF filtering."""
    import pandas as pd
    dfs = []
    for tree in trees:
        bed_file = get_position_bed_file(reference_genome, use_old, ancient_dna, tree)
        df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'name'])
        dfs.append(df)
    merged = pd.concat(dfs).drop_duplicates(subset='start').sort_values('start').reset_index(drop=True)
    merged['chr'] = 'chrY'  # main_vcf_split replaces this with the actual name from the VCF
    merged.to_csv(str(merged_bed), sep='\t', index=False, header=False)


def main_vcf_multi_tree(
        args: argparse.Namespace,
        base_out_folder: Path,
        trees: List[str]
):
    """Multi-tree VCF: filter once with union positions, extract per-tree, write combined output."""
    merged_bed = base_out_folder / "temp_vcf_union.bed"
    _write_merged_vcf_bed(merged_bed, trees, args.reference_genome, args.use_old, args.ancient_DNA)

    safe_create_dir(base_out_folder / "filtered_vcf_files", args.force)
    files = get_files_with_extension(args.vcffile, '.vcf.gz')

    if not args.reanalyze:
        with multiprocessing.Pool(processes=args.threads) as p:
            sample_vcf_files_nested = p.map(
                partial(main_vcf_split, merged_bed, base_out_folder, args), files)
        merged_bed.unlink(missing_ok=True)
        sample_vcf_files = [x for sublist in sample_vcf_files_nested if sublist for x in sublist]
    else:
        merged_bed.unlink(missing_ok=True)
        sample_vcf_files = files

    for tree in trees:
        path_markerfile = get_position_file(args.reference_genome, args.use_old, args.ancient_DNA, tree)
        with multiprocessing.Pool(processes=args.threads) as p:
            p.map(partial(run_vcf, path_markerfile, base_out_folder, args,
                          out_suffix=f'.{tree}', tree_name=tree),
                  sample_vcf_files)

    combined_out = base_out_folder / "hg_prediction_combined.hg"
    write_combined_prediction_table(base_out_folder, trees, combined_out,
                                    args.use_old, args.prediction_quality, args.threads,
                                    draw_hg=args.draw_haplogroups,
                                    collapsed_draw_mode=args.collapsed_draw_mode)


def main_vcf(
        args: argparse.Namespace,
        base_out_folder: Path
):
    position_bed_file = get_position_bed_file(args.reference_genome, args.use_old, args.ancient_DNA, args.tree)
    path_markerfile = get_position_file(args.reference_genome, args.use_old, args.ancient_DNA, args.tree)

    safe_create_dir(base_out_folder / "filtered_vcf_files", args.force)

    files = get_files_with_extension(args.vcffile, '.vcf.gz')

    if not args.reanalyze:
        with multiprocessing.Pool(processes=args.threads) as p:
            sample_vcf_files = p.map(partial(main_vcf_split, position_bed_file, base_out_folder, args), files)

        sample_vcf_files = [x for x in sample_vcf_files if x is not None]
        sample_vcf_files = [item for sublist in sample_vcf_files for item in sublist]

        with multiprocessing.Pool(processes=args.threads) as p:
            p.map(partial(run_vcf, path_markerfile, base_out_folder, args), sample_vcf_files)

    else:
        with multiprocessing.Pool(processes=args.threads) as p:
            p.map(partial(run_vcf, path_markerfile, base_out_folder, args), files)


def main():
    print("Erasmus MC Department of Genetic Identification\nYleaf: software tool for human Y-chromosomal "
          f"phylogenetic analysis and haplogroup inference v{__version__}")
    logo()

    args = get_arguments()
    out_folder = Path(args.output)
    safe_create_dir(out_folder, args.force, reuse_pileup=getattr(args, 'reuse_pileup', False))
    setup_logger(out_folder)

    LOG.info(f"Running Yleaf with command: {' '.join(sys.argv)}")

    # Normalize --tree: unwrap single-element list to string for backward compatibility
    if isinstance(args.tree, list) and len(args.tree) == 1:
        args.tree = args.tree[0]

    multi_tree = isinstance(args.tree, list)  # True when multiple trees were specified
    trees = args.tree if multi_tree else [args.tree]

    # validate tree/reference combinations
    for tree in trees:
        if tree == yleaf_constants.TREE_FTDNA:
            if args.use_old:
                LOG.error("--use_old is not compatible with --tree ftdna.")
                raise ValueError("--use_old is not compatible with --tree ftdna.")

    # make sure the reference genome is present before doing something else, if not present it is downloaded
    check_reference(args.reference_genome)

    if multi_tree:
        if args.fastq:
            main_fastq(args, out_folder, trees=trees)
        elif args.vcffile:
            main_vcf_multi_tree(args, out_folder, trees)
        elif args.plinkfile:
            main_plink_multi_tree(args, out_folder, trees)
        else:
            is_bam = args.bamfile is not None
            main_bam_cram_multi_tree(args, out_folder, is_bam, trees)
        LOG.info("Done!")
        return

    if args.fastq:
        main_fastq(args, out_folder)
    elif args.bamfile:
        main_bam_cram(args, out_folder, True)
    elif args.cramfile:
        main_bam_cram(args, out_folder, False)
    elif args.vcffile:
        main_vcf(args, out_folder)
    elif args.plinkfile:
        main_plink(args, out_folder)
    else:
        LOG.error("Please specify either a bam, cram, fastq, vcf, or plink file")
        raise ValueError("Please specify either a bam, cram, fastq, vcf, or plink file")
    hg_out = out_folder / PREDICTION_OUT_FILE_NAME
    predict_haplogroup(out_folder, hg_out, args.use_old, args.prediction_quality, args.threads, args.tree)
    _write_path_file(out_folder, args.tree, hg_out, ".out")
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
    parser.add_argument("-cr", "--cram_reference", required=False,
                        help="Reference genome for the CRAM file. Required when using CRAM files.",
                        metavar="PATH", type=check_file)
    parser.add_argument("-vcf", "--vcffile", required=False,
                        help="input VCF file (.vcf.gz)", metavar="PATH", type=check_file)
    parser.add_argument("-plink", "--plinkfile", required=False,
                        help="PLINK format SNP array input. Pass the .ped file (text format) or the "
                             ".bed file (binary format); companion files (.map/.bim/.fam) must be in "
                             "the same directory. Both single- and multi-sample files are supported. "
                             "Chromosome Y must be encoded as 24, Y, or chrY.",
                        metavar="PATH", type=check_file)
    parser.add_argument("-ra", "--reanalyze", required=False,
                        help="reanalyze (skip filtering and splitting) the vcf file", action="store_true")
    parser.add_argument("-rp", "--reuse_pileup", required=False,
                        help="reuse existing pileup file (.pu) if present, skipping samtools mpileup. "
                             "The .pu file is kept on disk when this flag is set.",
                        action="store_true")
    parser.add_argument("-force", "--force", action="store_true",
                        help="Delete files without asking")
    parser.add_argument("-rg", "--reference_genome",
                        help="The reference genome build to be used. If no reference is available "
                             "they will be downloaded. If you added references in your config.txt file these"
                             " will be used instead as reference or the location will be used to download the "
                             "reference if those files are missing or empty.",
                        choices=[yleaf_constants.HG19, yleaf_constants.HG38, yleaf_constants.T2T], required=True)
    parser.add_argument("-o", "--output", required=True,
                        help="Folder name containing outputs", metavar="STRING")
    parser.add_argument("-r", "--reads_treshold",
                        help="The minimum number of reads for each base. (default=10)",
                        type=int, required=False,
                        default=10, metavar="INT")
    parser.add_argument("-q", "--quality_thresh",
                        help="Minimum quality for each read, integer between 10 and 40. [10-40] (default=20)",
                        type=int, required=False, metavar="INT", default=20)
    parser.add_argument("-b", "--base_majority",
                        help="The minimum percentage of a base result for acceptance, integer between 50 and 99."
                             " [50-99] (default=90)",
                        type=int, required=False, metavar="INT", default=90)
    parser.add_argument("-t", "--threads", dest="threads",
                        help="The number of processes to use when running Yleaf.",
                        type=int, default=1, metavar="INT")
    parser.add_argument("-pq", "--prediction_quality", type=float, required=False, default=0.95, metavar="FLOAT",
                        help="The minimum quality of the prediction (QC-scores) for it to be accepted. [0-1] (default=0.95)")

    _tree_choices = [yleaf_constants.TREE_YFULL, yleaf_constants.TREE_YFULL_V10,
                     yleaf_constants.TREE_FTDNA, yleaf_constants.TREE_OPENYTREE,
                     yleaf_constants.TREE_ISOGG]
    parser.add_argument("-tree", "--tree",
                        help="One or more haplogroup trees to use for prediction. Accepted values: "
                             "yfull (YFull v14.01), yfull_v10 (YFull v10.01 legacy), "
                             "ftdna (FamilyTreeDNA, all references), openY (Open Y combined tree), "
                             "isogg (ISOGG legacy). When multiple trees are given a single combined "
                             "pileup is built and prediction is run per tree. (default=yfull)",
                        nargs='+',
                        choices=_tree_choices,
                        default=[yleaf_constants.TREE_YFULL])

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

    # arguments for ancient DNA samples
    parser.add_argument("-aDNA", "--ancient_DNA",
                        help="Add this flag if the sample is ancient DNA. This will ignore all G > A and C > T mutations.",
                        action="store_true")

    # arguments for private mutations
    parser.add_argument("-p", "--private_mutations",
                        help="Add this flag to search for private mutations. These are variations that are not"
                             " considered in the phylogenetic tree and thus not used for haplogroup prediction, "
                             "however can be informative and differentiate individuals within the same haplogroup "
                             "prediction.",
                        action="store_true")

    parser.add_argument("-maf", "--minor_allele_frequency", help="Maximum rate of minor allele for it to be considered"
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
    LOG.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    start_time = time.time()
    formatter = MyFormatter('%(levelname)s %(currentTime)s (%(passedTime)s s) - %(message)s',
                            starttime=start_time)
    handler.setFormatter(formatter)
    handler.setLevel(logging.INFO)
    LOG.addHandler(handler)

    file_handler = logging.FileHandler(filename=out_folder / "run.log")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    LOG.addHandler(file_handler)

    LOG.debug("Logger created")


def safe_create_dir(
        folder: Path,
        force: bool,
        reuse_pileup: bool = False
):
    """Create the given folder. If the folder is already present delete if the user agrees.
    When reuse_pileup is True, an existing folder is kept as-is (pileup reuse requires it)."""
    if folder.is_dir():
        if reuse_pileup:
            return  # keep existing folder so the pileup file can be found
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
            os.makedirs(folder, exist_ok=True)
        except OSError:
            print("Failed to create directory. Exiting...")
            raise


def main_fastq(
        args: argparse.Namespace,
        out_folder: Path,
        trees: List[str] = None
):
    files = get_files_with_extension(args.fastq, '.fastq')
    files += get_files_with_extension(args.fastq, '.fastq.gz')
    reference = get_reference_path(args.reference_genome, True)
    bam_folder = out_folder / yleaf_constants.FASTQ_BAM_FILE_FOLDER
    try:
        os.mkdir(bam_folder)
    except IOError:
        pass
    LOG.info("Creating bam files from fastq files...")

    # For paired end fastq files (with _R1 and _R2 in the file name) and gzipped fastq files with .gz extension align the pairs together
    for fastq_file in files:
        if "_R1" in str(fastq_file) and Path(str(fastq_file).replace("_R1", "_R2")).exists():
            fastq_file2 = Path(str(fastq_file).replace("_R1", "_R2"))
            LOG.info(f"Starting with running for {fastq_file} and {fastq_file2}")
            sam_file = bam_folder / "temp_fastq_sam.sam"
            fastq_cmd = f"minimap2 -ax sr -k14 -w7 -t {args.threads} {reference} {fastq_file} {fastq_file2} > {sam_file}"
            call_command(fastq_cmd)
            bam_file = bam_folder / (fastq_file.name.rsplit("_R1", 1)[0] + ".bam")
            cmd = "samtools view -@ {} -bS {} | samtools sort -@ {} -m 2G -o {}".format(args.threads, sam_file,
                                                                                        args.threads, bam_file)
            call_command(cmd)
            cmd = "samtools index -@ {} {}".format(args.threads, bam_file)
            call_command(cmd)
            os.remove(sam_file)
        elif "_R1" not in str(fastq_file) and "_R2" not in str(fastq_file) and ".gz" in str(fastq_file):
            LOG.info(f"Starting with running for {fastq_file}")
            sam_file = bam_folder / "temp_fastq_sam.sam"

            fastq_cmd = f"minimap2 -ax sr -k14 -w7 -t {args.threads} {reference} {fastq_file} > {sam_file}"
            call_command(fastq_cmd)
            bam_file = bam_folder / (fastq_file.name.rsplit(".", 1)[0] + ".bam")
            cmd = "samtools view -@ {} -bS {} | samtools sort -@ {} -m 2G -o {}".format(args.threads, sam_file,
                                                                                        args.threads, bam_file)
            call_command(cmd)
            cmd = "samtools index -@ {} {}".format(args.threads, bam_file)
            call_command(cmd)
            os.remove(sam_file)
    args.bamfile = bam_folder
    if trees is not None:
        main_bam_cram_multi_tree(args, out_folder, True, trees)
    else:
        main_bam_cram(args, out_folder, True)


_CRAM_CHRY_MD5 = {
    "ce3e31103314a704255f3cd90369ecce": yleaf_constants.HG38,
    "17ceeafc3a372967496cd94f2572c341": yleaf_constants.HG19,
    "360f2d39c52a70539bac6b12ff76c123": yleaf_constants.T2T,
}


def resolve_cram_reference(args: argparse.Namespace, cram_files: List[Path]) -> None:
    """Auto-detect reference genome from CRAM header chrY M5 checksum and set args.cram_reference."""
    if args.cram_reference is not None:
        return
    if not cram_files:
        raise ValueError("No CRAM files found.")

    result = subprocess.run(
        f"samtools view -H {cram_files[0]}",
        shell=True, capture_output=True, text=True
    )
    detected_build = None
    for line in result.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
        if fields.get("SN") in ("chrY", "Y") and "M5" in fields:
            detected_build = _CRAM_CHRY_MD5.get(fields["M5"])
            if detected_build:
                break

    if detected_build is None:
        raise ValueError(
            "Cannot auto-detect reference genome from CRAM header (no recognised chrY M5 checksum). "
            "Please specify the reference FASTA with -cr/--cram_reference."
        )

    LOG.info(f"Auto-detected CRAM reference genome from chrY M5: {detected_build}")

    ref_path = {
        yleaf_constants.HG38: yleaf_constants.HG38_FULL_GENOME,
        yleaf_constants.HG19: yleaf_constants.HG19_FULL_GENOME,
        yleaf_constants.T2T:  yleaf_constants.T2T_FULL_GENOME,
    }[detected_build]

    if not ref_path.exists():
        LOG.info(f"Reference genome not found at {ref_path}. Downloading (~3 GB), please wait...")
        from yleaf import download_reference
        download_reference.install_genome_files(detected_build)

    if args.reference_genome != detected_build:
        LOG.warning(
            f"Auto-detected reference build ({detected_build}) differs from -rg "
            f"({args.reference_genome}); using detected build."
        )
        args.reference_genome = detected_build
    args.cram_reference = ref_path


def main_bam_cram(
        args: argparse.Namespace,
        base_out_folder: Path,
        is_bam: bool
):
    if args.bamfile is not None:
        files = get_files_with_extension(args.bamfile, '.bam')
    elif args.cramfile is not None:
        files = get_files_with_extension(args.cramfile, '.cram')
        resolve_cram_reference(args, files)
    else:
        print("Please specify either (a) bam or a cram file(s)")
        return

    with multiprocessing.Pool(processes=args.threads) as p:
        p.map(partial(run_bam_cram, args, base_out_folder, is_bam), files)


def main_bam_cram_multi_tree(
        args: argparse.Namespace,
        base_out_folder: Path,
        is_bam: bool,
        trees: List[str]
):
    """Run multi-tree analysis: one merged pileup per sample, prediction per tree, combined output."""
    if args.bamfile is not None:
        files = get_files_with_extension(args.bamfile, '.bam')
    elif args.cramfile is not None:
        files = get_files_with_extension(args.cramfile, '.cram')
        resolve_cram_reference(args, files)
    else:
        return

    with multiprocessing.Pool(processes=args.threads) as p:
        p.map(partial(run_bam_cram_multi_tree, args, base_out_folder, is_bam, trees), files)

    # Flush all pending NFS/OS writes from worker processes before reading output files
    os.sync()

    # Verify all expected per-tree .out files are present and readable.
    # On CIFS mounts, reading content (not just stat) forces the page cache to be
    # populated before prediction workers fork, preventing stale-read issues.
    missing = []
    for f in files:
        sample_name = f.name.rsplit(".", 1)[0]
        sample_dir = base_out_folder / sample_name
        for tree in trees:
            out_file = sample_dir / f"{sample_name}.{tree}.out"
            try:
                with open(out_file) as fh:
                    header = fh.readline()
                if not header.startswith("chr"):
                    missing.append(str(out_file) + " (bad header)")
            except FileNotFoundError:
                missing.append(str(out_file) + " (missing)")
    if missing:
        LOG.error(f"{len(missing)} expected .out file(s) missing or unreadable after extraction — "
                  f"prediction results may be incomplete:\n" + "\n".join(missing))

    # Write combined prediction table (and optionally draw per-tree haplogroup images)
    combined_out = base_out_folder / "hg_prediction_combined.hg"
    write_combined_prediction_table(base_out_folder, trees, combined_out,
                                    args.use_old, args.prediction_quality, args.threads,
                                    draw_hg=args.draw_haplogroups,
                                    collapsed_draw_mode=args.collapsed_draw_mode)


def run_bam_cram_multi_tree(
        args: argparse.Namespace,
        base_out_folder: Path,
        is_bam: bool,
        trees: List[str],
        input_file: Path
):
    """Per-sample multi-tree: one pileup covering all trees, extract haplogroups per tree."""
    try:
        _run_bam_cram_multi_tree_inner(args, base_out_folder, is_bam, trees, input_file)
    except BaseException as e:
        LOG.warning(f"Skipping {input_file.name}: {e}")


def _run_bam_cram_multi_tree_inner(
        args: argparse.Namespace,
        base_out_folder: Path,
        is_bam: bool,
        trees: List[str],
        input_file: Path
):
    LOG.info(f"Starting multi-tree run for {input_file}")
    output_dir = base_out_folder / input_file.name.rsplit(".", 1)[0]
    safe_create_dir(output_dir, args.force, reuse_pileup=getattr(args, 'reuse_pileup', False))

    if is_bam:
        if not any([Path(str(input_file) + ".bai").exists(),
                    Path(str(input_file).rstrip(".bam") + '.bai').exists()]):
            call_command(f"samtools index -@ {args.threads} {input_file}")
    else:
        if not any([Path(str(input_file) + ".crai").exists(),
                    Path(str(input_file).rstrip(".cram") + '.crai').exists()]):
            call_command(f"samtools index -@ {args.threads} {input_file}")

    header, mapped, unmapped = chromosome_table(input_file, output_dir, output_dir.name)
    general_info_list = [f"Total of mapped reads: {mapped}", f"Total of unmapped reads: {unmapped}"]

    # Build merged BED from all trees
    pileupfile = output_dir / "temp_haplogroup_pileup.pu"
    reuse_pileup = getattr(args, 'reuse_pileup', False)

    if reuse_pileup and pileupfile.exists():
        LOG.info(f"Reusing existing pileup file: {pileupfile}")
    else:
        merged_bed = output_dir / "temp_position_bed_combined.bed"
        _write_merged_bed(merged_bed, trees, args.reference_genome, args.use_old, args.ancient_DNA, header)
        reference = args.cram_reference if not is_bam else None
        execute_mpileup(merged_bed, input_file, pileupfile, args.quality_thresh, reference)
        os.remove(merged_bed)

    # Extract haplogroups per tree, writing a per-tree .info file for each
    for tree in trees:
        position_file = get_position_file(args.reference_genome, args.use_old, args.ancient_DNA, tree)
        out_suffix = f".{tree}"
        fmf_output = output_dir / (output_dir.name + out_suffix + ".fmf")
        outputfile = output_dir / (output_dir.name + out_suffix + ".out")
        tree_info_list = list(general_info_list)  # fresh copy per tree
        extract_haplogroups(position_file, args.reads_treshold, args.base_majority,
                            pileupfile, fmf_output, outputfile, is_bam, args.use_old,
                            tree_info_list, tree)
        write_info_file(output_dir, tree_info_list, suffix=f".{tree}.info")

    # TODO (production): delete pileup when no longer needed: os.remove(pileupfile)

    LOG.debug(f"Finished multi-tree extraction for {input_file}")


def _write_merged_bed(
        merged_bed: Path,
        trees: List[str],
        reference_genome: str,
        use_old: bool,
        ancient_dna: bool,
        header: str
):
    """Merge BED files from multiple trees, deduplicate positions, write sorted output."""
    import pandas as pd
    dfs = []
    for tree in trees:
        pos_file = get_position_file(reference_genome, use_old, ancient_dna, tree)
        df = pd.read_csv(pos_file, sep='\t', header=None, usecols=[0, 3])
        df.columns = ['chr', 'pos']
        df['chr'] = header
        dfs.append(df)
    merged = pd.concat(dfs).drop_duplicates(subset='pos').sort_values('pos').reset_index(drop=True)
    merged.to_csv(str(merged_bed), sep='\t', index=False, header=False)


def write_combined_prediction_table(
        base_out_folder: Path,
        trees: List[str],
        combined_out: Path,
        use_old: bool,
        prediction_quality: float,
        threads: int,
        draw_hg: bool = False,
        collapsed_draw_mode: bool = False
):
    """Run predict_haplogroup per tree and write combined TSV with one row per sample."""
    # Collect per-tree predictions
    tree_results = {}  # tree -> {sample_name -> (hg, hg_marker, valid_markers, qc)}
    for tree in trees:
        tree_hg_out = base_out_folder / f"hg_prediction_{tree}.hg"
        _predict_for_tree(base_out_folder, tree, tree_hg_out, use_old, prediction_quality, threads)
        tree_results[tree] = _read_hg_file(tree_hg_out)
        if draw_hg:
            draw_haplogroups(
                tree_hg_out, collapsed_draw_mode,
                outfile=base_out_folder / f"hg_tree_image_{tree}",
                tree_file=get_tree_path(tree)
            )

    # Collect all sample names
    sample_names = set()
    for results in tree_results.values():
        sample_names.update(results.keys())

    with open(combined_out, 'w') as f:
        header_cols = ['Sample_name']
        for tree in trees:
            header_cols += [f'{tree}_Hg', f'{tree}_Hg_marker', f'{tree}_Valid_markers', f'{tree}_QC']
        f.write('\t'.join(header_cols) + '\n')
        for sample in sorted(sample_names):
            row = [sample]
            for tree in trees:
                result = tree_results[tree].get(sample, ('NA', 'NA', 'NA', 'NA'))
                row += list(result)
            f.write('\t'.join(str(x) for x in row) + '\n')

    LOG.info(f"Combined prediction table written to {combined_out}")


def _predict_for_tree(
        base_out_folder: Path,
        tree: str,
        hg_out: Path,
        use_old: bool,
        prediction_quality: float,
        threads: int
):
    """Run haplogroup prediction for one tree, reading tree-specific .{tree}.out files."""
    from yleaf import predict_haplogroup
    namespace = argparse.Namespace(
        input=base_out_folder,
        outfile=hg_out,
        minimum_score=prediction_quality,
        threads=threads,
        tree_file=get_tree_path(tree),
        out_file_suffix=f".{tree}.out",
        info_file_suffix=f".{tree}.info",
    )
    predict_haplogroup.main(namespace)
    _write_path_file(base_out_folder, tree, hg_out, f".{tree}.out")


def _write_path_file(
        base_out_folder: Path,
        tree_name: str,
        hg_out_file: Path,
        out_file_suffix: str,
):
    """Write per-sample haplogroup_path_{tree}.json inside each sample's subfolder."""
    import json as _json
    if not hg_out_file.exists():
        return

    # Read predicted haplogroup for every sample from the .hg file
    sample_hg: dict = {}
    with open(hg_out_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2 and parts[1] not in ('NA', ''):
                sample_hg[parts[0]] = parts[1]

    if not sample_hg:
        return

    tree = Tree(get_tree_path(tree_name))

    for sample_name, predicted_hg in sample_hg.items():
        # Strip subclade annotation e.g. "E-Z15929*(xE-Y25504)" → "E-Z15929"
        predicted_hg_base = predicted_hg.split('*')[0].split('(')[0]
        node = tree.get(predicted_hg_base)
        if not node:
            continue

        # Read marker states from this sample's own .out file
        sample_folder = base_out_folder / sample_name
        out_path = sample_folder / f"{sample_name}{out_file_suffix}"
        hg_counts: dict = {}
        if out_path.exists():
            with open(out_path) as f:
                next(f)
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 11:
                        hg, marker, state = parts[3], parts[2], parts[10]
                        if hg not in hg_counts:
                            hg_counts[hg] = {'d': 0, 'a': 0, 'rep_d': None, 'rep_a': None}
                        if state == 'D':
                            hg_counts[hg]['d'] += 1
                            if hg_counts[hg]['rep_d'] is None:
                                hg_counts[hg]['rep_d'] = marker
                        elif state == 'A':
                            hg_counts[hg]['a'] += 1
                            if hg_counts[hg]['rep_a'] is None:
                                hg_counts[hg]['rep_a'] = marker

        marker_state: dict = {}
        for hg, v in hg_counts.items():
            majority = 'D' if v['d'] >= v['a'] else 'A'
            rep = v['rep_d'] if majority == 'D' else v['rep_a']
            marker_state[hg] = {'marker': rep or v['rep_d'] or v['rep_a'], 'state': majority}

        path = []
        n = node
        while n:
            entry = marker_state.get(n.name)
            path.append({
                'haplogroup': n.name,
                'marker': entry['marker'] if entry else None,
                'state': entry['state'] if entry else 'NC',
            })
            n = n.parent
        path.reverse()

        out = sample_folder / f'haplogroup_path_{tree_name}.json'
        with open(out, 'w') as f:
            _json.dump(path, f)


def _read_hg_file(hg_file: Path) -> dict:
    """Read a hg_prediction.hg file, return {sample_name: (Hg, Hg_marker, Valid_markers, QC)}."""
    results = {}
    if not hg_file.exists():
        return results
    with open(hg_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                results[parts[0]] = (parts[1], parts[2], parts[4], parts[5])
    return results


def run_bam_cram(
        args: argparse.Namespace,
        base_out_folder: Path,
        is_bam: bool,
        input_file: Path
):
    try:
        LOG.info(f"Starting with running for {input_file}")
        output_dir = base_out_folder / input_file.name.rsplit(".", 1)[0]
        safe_create_dir(output_dir, args.force, reuse_pileup=getattr(args, 'reuse_pileup', False))
        general_info_list = samtools(output_dir, input_file, is_bam, args)
        write_info_file(output_dir, general_info_list)
        if args.private_mutations:
            find_private_mutations(output_dir, input_file, args, is_bam)
        LOG.debug(f"Finished running for {input_file.name}")
        print()
    except BaseException as e:
        LOG.warning(f"Skipping {input_file.name}: {e}")


def call_command(
        command_str: str,
        stdout_location: TextIO = None
):
    """Call a command on the command line and make sure to exit if it fails."""
    LOG.debug(f"Started running the following command: {command_str}")
    if stdout_location is None:
        process = subprocess.Popen(command_str, stderr=subprocess.PIPE, shell=True)
    else:
        process = subprocess.Popen(command_str, stderr=subprocess.PIPE, stdout=stdout_location, shell=True)
    # blocking call
    stdout, stderr = process.communicate()
    # will only fail if returncode is not 0
    if process.returncode != 0:
        LOG.error(f"Call: '{command_str}' failed. Reason given: '{stderr.decode('utf-8')}'")
        raise SystemExit("Failed command execution")
    LOG.debug(f"Finished running the command {command_str}")


def get_files_with_extension(
        path: Union[str, Path],
        ext: str
) -> List[Path]:
    """Get all files with a certain extension from a path. The path can be a file or a dir."""
    filtered_files = []
    path = Path(path)  # to be sure
    if path.is_dir():
        for file in path.iterdir():
            if str(file)[-len(ext):] == ext:
                filtered_files.append(file)
        return filtered_files
    else:
        return [path]


def check_reference(
        requested_version: str,
):
    reference_file = get_reference_path(requested_version, True)
    if not reference_file.exists() or os.path.getsize(reference_file) < 100:
        LOG.info(f"No reference genome version was found. Downloading the {requested_version} reference genome. This "
                 f"should be a one time thing.")
        try:
            download_reference.main(requested_version)
        except Exception as e:
            LOG.error(
                f"Failed to download the {requested_version} reference genome: {e}\n"
                f"Please check your internet connection and ensure that the data directory "
                f"({reference_file.parent}) is writable, then re-run Yleaf. "
                f"Alternatively, place the reference FASTA manually and set its path in "
                f"yleaf/config.txt."
            )
            sys.exit(1)
        LOG.info("Finished downloading the reference genome.")


def get_reference_path(
        requested_version: str,
        is_full: bool
) -> Union[Path, None]:
    if is_full:
        if requested_version == yleaf_constants.HG19:
            reference_file = yleaf_constants.HG19_FULL_GENOME
        elif requested_version == yleaf_constants.T2T:
            reference_file = yleaf_constants.T2T_FULL_GENOME
        else:
            reference_file = yleaf_constants.HG38_FULL_GENOME
    else:
        if requested_version == yleaf_constants.HG19:
            reference_file = yleaf_constants.HG19_Y_CHROMOSOME
        elif requested_version == yleaf_constants.T2T:
            reference_file = yleaf_constants.T2T_Y_CHROMOSOME
        else:
            reference_file = yleaf_constants.HG38_Y_CHROMOSOME
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
    reference = args.cram_reference if not is_bam_pathfile else None

    if is_bam_pathfile:
        if not any([Path(str(path_file) + ".bai").exists(), Path(str(path_file).rstrip(".bam") + '.bai').exists()]):
            cmd = "samtools index -@{} {}".format(args.threads, path_file)
            call_command(cmd)
    else:
        if not any([Path(str(path_file) + ".crai").exists(), Path(str(path_file).rstrip(".cram") + '.crai').exists()]):
            cmd = "samtools index -@{} {}".format(args.threads, path_file)
            call_command(cmd)
    header, mapped, unmapped = chromosome_table(path_file, output_folder, output_folder.name)

    position_file = get_position_file(args.reference_genome, args.use_old, args.ancient_DNA, args.tree)

    reuse_pileup = getattr(args, 'reuse_pileup', False)
    if reuse_pileup and pileupfile.exists():
        LOG.info(f"Reusing existing pileup file: {pileupfile}")
        bed = None
    else:
        bed = output_folder / "temp_position_bed.bed"
        write_bed_file(bed, position_file, header)
        execute_mpileup(bed, path_file, pileupfile, args.quality_thresh, reference)

    general_info_list = ["Total of mapped reads: " + str(mapped), "Total of unmapped reads: " + str(unmapped)]

    extract_haplogroups(position_file, args.reads_treshold, args.base_majority,
                        pileupfile, fmf_output, outputfile, is_bam_pathfile, args.use_old, general_info_list,
                        args.tree)

    # TODO (production): delete pileup when no longer needed: os.remove(pileupfile)
    if bed is not None:
        os.remove(bed)

    LOG.debug("Finished extracting haplogroups")
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


def get_tree_path(tree: str) -> Path:
    if tree == yleaf_constants.TREE_FTDNA:
        return yleaf_constants.HG_PREDICTION_FOLDER / yleaf_constants.FTDNA_TREE_FILE
    if tree == yleaf_constants.TREE_OPENYTREE:
        return yleaf_constants.HG_PREDICTION_FOLDER / yleaf_constants.OPENYTREE_TREE_FILE
    if tree == yleaf_constants.TREE_ISOGG:
        return yleaf_constants.HG_PREDICTION_FOLDER / yleaf_constants.ISOGG_TREE_FILE
    if tree == yleaf_constants.TREE_YFULL_V10:
        return yleaf_constants.HG_PREDICTION_FOLDER / yleaf_constants.YFULL_V10_TREE_FILE
    return yleaf_constants.HG_PREDICTION_FOLDER / yleaf_constants.TREE_FILE


def get_position_file(
        reference_name: str,
        use_old: bool,
        ancient_DNA: bool,
        tree: str = yleaf_constants.TREE_YFULL,
) -> Path:
    if tree == yleaf_constants.TREE_FTDNA:
        if ancient_DNA:
            fname = yleaf_constants.FTDNA_POSITION_ANCIENT_FILE.format(ref=reference_name)
        else:
            fname = yleaf_constants.FTDNA_POSITION_FILE.format(ref=reference_name)
        return yleaf_constants.DATA_FOLDER / reference_name / fname
    if tree == yleaf_constants.TREE_OPENYTREE:
        ref_tag = {"hg38": "hg38", "hg19": "hg19", "t2t": "t2t"}.get(reference_name, reference_name)
        fname = yleaf_constants.OPENYTREE_POSITION_FILE.format(ref=ref_tag)
        return yleaf_constants.DATA_FOLDER / reference_name / fname
    if tree == yleaf_constants.TREE_ISOGG:
        ref_tag = {"hg38": "hg38", "hg19": "hg19", "t2t": "t2t"}.get(reference_name, reference_name)
        fname = yleaf_constants.ISOGG_POSITION_FILE.format(ref=ref_tag)
        return yleaf_constants.DATA_FOLDER / reference_name / fname
    if tree == yleaf_constants.TREE_YFULL_V10:
        ref_tag = {"hg38": "hg38", "hg19": "hg19", "t2t": "t2t"}.get(reference_name, reference_name)
        fname = yleaf_constants.YFULL_V10_POSITION_FILE.format(ref=ref_tag)
        return yleaf_constants.DATA_FOLDER / reference_name / fname
    if use_old:
        if ancient_DNA:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.OLD_POSITION_ANCIENT_FILE
        else:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.OLD_POSITION_FILE
    else:
        if ancient_DNA:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.NEW_POSITION_ANCIENT_FILE
        else:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.NEW_POSITION_FILE
    return position_file


def get_position_bed_file(
        reference_name: str,
        use_old: bool,
        ancient_DNA: bool,
        tree: str = yleaf_constants.TREE_YFULL,
) -> Path:
    if tree == yleaf_constants.TREE_FTDNA:
        if ancient_DNA:
            fname = yleaf_constants.FTDNA_POSITION_ANCIENT_BED_FILE.format(ref=reference_name)
        else:
            fname = yleaf_constants.FTDNA_POSITION_BED_FILE.format(ref=reference_name)
        return yleaf_constants.DATA_FOLDER / reference_name / fname
    if tree == yleaf_constants.TREE_OPENYTREE:
        ref_tag = {"hg38": "hg38", "hg19": "hg19", "t2t": "t2t"}.get(reference_name, reference_name)
        fname = yleaf_constants.OPENYTREE_POSITION_BED_FILE.format(ref=ref_tag)
        return yleaf_constants.DATA_FOLDER / reference_name / fname
    if tree == yleaf_constants.TREE_ISOGG:
        ref_tag = {"hg38": "hg38", "hg19": "hg19", "t2t": "t2t"}.get(reference_name, reference_name)
        fname = yleaf_constants.ISOGG_POSITION_BED_FILE.format(ref=ref_tag)
        return yleaf_constants.DATA_FOLDER / reference_name / fname
    if tree == yleaf_constants.TREE_YFULL_V10:
        ref_tag = {"hg38": "hg38", "hg19": "hg19", "t2t": "t2t"}.get(reference_name, reference_name)
        fname = yleaf_constants.YFULL_V10_POSITION_BED_FILE.format(ref=ref_tag)
        return yleaf_constants.DATA_FOLDER / reference_name / fname
    if use_old:
        if ancient_DNA:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.OLD_POSITION_ANCIENT_BED_FILE
        else:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.OLD_POSITION_BED_FILE
    else:
        if ancient_DNA:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.NEW_POSITION_ANCIENT_BED_FILE
        else:
            position_file = yleaf_constants.DATA_FOLDER / reference_name / yleaf_constants.NEW_POSITION_BED_FILE
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
    cmd += f" -ABQ{quality_thresh}q1 {str(bam_file)} > {str(pileupfile)}"
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
        general_info_list: List[str],
        tree: str = yleaf_constants.TREE_YFULL,
):
    LOG.debug("Starting with extracting haplogroups...")
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

    tree = Tree(get_tree_path(tree))
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
        if char in {"-", "+"}:
            # Indel notation: +<n><bases> or -<n><bases>. Count the marker but
            # skip the following bases so they are not counted as real reads.
            # This check must come before the fastadict check because "-" and "+"
            # are also keys in fastadict (issue #31).
            fastadict[char] += 1
            index += 1
            digit, index = find_digit(sequence, index)
            index += digit
        elif char == "^":
            index += 2
        elif char in fastadict:
            fastadict[char] += 1
            index += 1
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
        general_info_list: List[str],
        suffix: str = ".info"
):
    try:
        with open(folder / (folder.name + suffix), "a") as f:
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
    LOG.debug("Starting with extracting private mutations...")
    snp_reference_file = get_reference_path(args.reference_genome, False)
    snp_database_file = yleaf_constants.DATA_FOLDER / args.reference_genome / yleaf_constants.SNP_DATA_FILE

    # run mpileup
    pileup_file = output_folder / "temp_private_mutation_pileup.pu"
    execute_mpileup(None, path_file, pileup_file, args.quality_thresh, snp_reference_file if not is_bam else None)

    LOG.debug("Loading reference files")

    position_file = get_position_file(args.reference_genome, args.use_old, args.ancient_DNA, args.tree)
    filter_positions = load_filter_data(position_file)
    known_snps = load_snp_database_file(snp_database_file, args.minor_allele_frequency)
    ychrom_reference = load_reference_file(snp_reference_file)

    LOG.debug("Finding private mutations...")
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
    with open(output_folder / f"{output_folder.name}.pmu", "w") as f:
        f.write(f"chrom\tposition\trn_no\tmutation\treference\tdetected\treads\tcalled_percentage\t"
                f"minor allele frequency\n")
        f.write(''.join(confirmed_private_mutations))
        f.write(''.join(private_mutations))
    LOG.debug("Finished extracting private mutations")


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
        minor_allele_frequency: float
) -> Union[Dict[str, List[Dict[str, str]]], None]:
    global CACHED_SNP_DATABASE

    if path is not None:
        if CACHED_SNP_DATABASE is None:
            CACHED_SNP_DATABASE = defaultdict(list)
            with open(path) as f:
                f.readline()
                for line in f:
                    rs, position, major_allele, minor_allele, frequency = line.strip().split(",")
                    if float(frequency) > minor_allele_frequency:
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
        use_old: bool,
        prediction_quality: float,
        threads: int,
        tree: str = yleaf_constants.TREE_YFULL,
):
    if use_old:
        script = yleaf_constants.SRC_FOLDER / "old_predict_haplogroup.py"
        cmd = "python {} -i {} -o {}".format(script, path_file, output)
        call_command(cmd)
    else:
        from yleaf import predict_haplogroup
        namespace = argparse.Namespace(input=path_file, outfile=output,
                                       minimum_score=prediction_quality, threads=threads,
                                       tree_file=get_tree_path(tree))
        predict_haplogroup.main(namespace)


def draw_haplogroups(
        haplogroup_file: Path,
        collapsed_draw_mode: bool,
        outfile: Path = None,
        tree_file: Path = None
):
    # make sure that it is only imported if requested by user
    from yleaf import haplogroup_tree_image
    if outfile is None:
        outfile = haplogroup_file.parent / HAPLOGROUP_IMAGE_FILE_NAME
    namespace = argparse.Namespace(input=haplogroup_file, collapse_mode=collapsed_draw_mode,
                                   outfile=outfile, tree_file=tree_file)
    haplogroup_tree_image.main(namespace)


if __name__ == "__main__":
    main()
