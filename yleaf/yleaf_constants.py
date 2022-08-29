from pathlib import Path


SRC_FOLDER: Path = Path(__file__).absolute().parent
DATA_FOLDER: Path = SRC_FOLDER / "data"
HG_PREDICTION_FOLDER: Path = DATA_FOLDER / "hg_prediction_tables"

HG19: str = "hg19"
HG38: str = "hg38"
HG19_FOLDER: Path = DATA_FOLDER / HG19
HG38_FOLDER: Path = DATA_FOLDER / HG38

FULL_REF_FILE: str = "full_reference.fa"
Y_REF_FILE: str = "chrY.fa"
SNP_DATA_FILE: str = "snp_data.csv"
NEW_POSITION_FILE: str = "new_positions.txt"
OLD_POSITION_FILE: str = "old_positions.txt"
TREE_FILE: str = "tree.json"

