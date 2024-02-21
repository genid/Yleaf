from pathlib import Path


SRC_FOLDER: Path = Path(__file__).absolute().parent
DATA_FOLDER: Path = SRC_FOLDER / "data"
CONFIG_PATH: Path = SRC_FOLDER / "config.txt"
CHR_NAMING_CONVENTION_FILE: Path = SRC_FOLDER / "chr_naming_convention.txt"
HG_PREDICTION_FOLDER: Path = DATA_FOLDER / "hg_prediction_tables"

HG19: str = "hg19"
HG38: str = "hg38"
__HG19_FOLDER: Path = DATA_FOLDER / HG19
__HG38_FOLDER: Path = DATA_FOLDER / HG38

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

TREE_FILE: str = "tree.json"

HG19_FULL_GENOME: Path = __HG19_FOLDER / FULL_REF_FILE
HG19_Y_CHROMOSOME: Path = __HG19_FOLDER / Y_REF_FILE
HG38_FULL_GENOME: Path = __HG38_FOLDER / FULL_REF_FILE
HG38_Y_CHROMOSOME: Path = __HG38_FOLDER / Y_REF_FILE


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


FASTQ_BAM_FILE_FOLDER: str = "bam_files"