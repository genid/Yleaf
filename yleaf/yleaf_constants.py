from pathlib import Path


SRC_FOLDER: Path = Path(__file__).absolute().parent
DATA_FOLDER: Path = SRC_FOLDER / "data"
CONFIG_PATH: Path = SRC_FOLDER / "config.txt"
CHR_NAMING_CONVENTION_FILE: Path = SRC_FOLDER / "chr_naming_convention.txt"
HG_PREDICTION_FOLDER: Path = DATA_FOLDER / "hg_prediction_tables"

HG19: str = "hg19"
HG38: str = "hg38"
T2T: str = "t2t"
__HG19_FOLDER: Path = DATA_FOLDER / HG19
__HG38_FOLDER: Path = DATA_FOLDER / HG38
__T2T_FOLDER: Path = DATA_FOLDER / T2T

FULL_REF_FILE: str = "full_reference.fa"
Y_REF_FILE: str = "chrY.fa"
SNP_DATA_FILE: str = "snp_data.csv"
NEW_POSITION_FILE: str = "new_positions.txt"
OLD_POSITION_FILE: str = "old_positions.txt"
NEW_POSITION_BED_FILE: str = "new_positions.bed"
OLD_POSITION_BED_FILE: str = "old_positions.bed"
NEW_POSITION_ANCIENT_FILE: str = "new_positions_ancient.txt"
OLD_POSITION_ANCIENT_FILE: str = "old_positions_ancient.txt"
NEW_POSITION_ANCIENT_BED_FILE: str = "new_positions_ancient.bed"
OLD_POSITION_ANCIENT_BED_FILE: str = "old_positions_ancient.bed"

TREE_YFULL: str = "yfull"
TREE_YFULL_V10: str = "yfull_v10"
TREE_FTDNA: str = "ftdna"
TREE_OPENYTREE: str = "openY"
TREE_ISOGG: str = "isogg"

TREE_FILE: str = "tree.json"
FTDNA_TREE_FILE: str = "ftdna_tree.json"
FTDNA_POSITION_FILE: str = "ftdna_positions_{ref}.txt"
FTDNA_POSITION_BED_FILE: str = "ftdna_positions_{ref}.bed"
OPENYTREE_TREE_FILE: str = "openY_tree.json"
OPENYTREE_POSITION_FILE: str = "openY_positions_{ref}.txt"
OPENYTREE_POSITION_BED_FILE: str = "openY_positions_{ref}.bed"
ISOGG_TREE_FILE: str = "isogg_tree.json"
ISOGG_POSITION_FILE: str = "isogg_positions_{ref}.txt"
ISOGG_POSITION_BED_FILE: str = "isogg_positions_{ref}.bed"
YFULL_V10_TREE_FILE: str = "yfull_v10_tree.json"
YFULL_V10_POSITION_FILE: str = "yfull_v10_positions_{ref}.txt"
YFULL_V10_POSITION_BED_FILE: str = "yfull_v10_positions_{ref}.bed"

HG19_FULL_GENOME: Path = __HG19_FOLDER / FULL_REF_FILE
HG19_Y_CHROMOSOME: Path = __HG19_FOLDER / Y_REF_FILE
HG38_FULL_GENOME: Path = __HG38_FOLDER / FULL_REF_FILE
HG38_Y_CHROMOSOME: Path = __HG38_FOLDER / Y_REF_FILE
T2T_FULL_GENOME: Path = __T2T_FOLDER / FULL_REF_FILE
T2T_Y_CHROMOSOME: Path = __T2T_FOLDER / Y_REF_FILE


def get_path(
    name_: str,
    value_: str
) -> Path:
    path = Path(value_)
    if not path.exists():
        if not path.parent.exists():
            raise ValueError(f"Cant find provided config path ({path}) for: '{name_}'. Try to define an absolute"
                             "path!")
        else:
            # create an empty file at the location
            open(path, "w").close()
    if path.suffix not in [".fa", ".fasta", ".fna"]:
        raise ValueError("Please provide a fasta file. File ending in .fa, .fasta or .fna")
    return path


# read the config to see if the paths are defined
with open(CONFIG_PATH) as f:
    for line in f:
        name, value = line.strip().split('=')
        value = value.strip()
        name = name.strip()
        if value == "":
            continue
        if name == "full hg19 genome fasta location":
            HG19_FULL_GENOME = get_path(name, value)
        elif name == "full hg38 genome fasta location":
            HG38_FULL_GENOME = get_path(name, value)
        elif name == "hg19 chromosome Y fasta location":
            HG19_Y_CHROMOSOME = get_path(name, value)
        elif name == "hg38 chromosome Y fasta location":
            HG38_Y_CHROMOSOME = get_path(name, value)
        elif name == "full t2t genome fasta location":
            T2T_FULL_GENOME = get_path(name, value)
        elif name == "t2t chromosome Y fasta location":
            T2T_Y_CHROMOSOME = get_path(name, value)


FASTQ_BAM_FILE_FOLDER: str = "bam_files"