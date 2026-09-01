<img src="yleaf_logo.png" width="500" alt="Yleaf logo">

# Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data

#### Arwin Ralf, Diego Montiel Gonzalez, Kaiyin Zhong and Manfred Kayser

### Department of Genetic Identification 
#### Erasmus MC University Medical Centre Rotterdam, The Netherlands

## Dashboard (GUI) — recommended for most users

The easiest way to use Yleaf is through the graphical dashboard. Download the installer for your platform from the [releases page](https://github.com/genid/Yleaf/releases):

**Windows**
- Download `Yleaf.4.1_4.1.5_x64-setup.exe` (or `Yleaf.4.1_4.1.5_x64_en-US.msi`) and run the installer.
- Launch **Yleaf 4.1** from the Start menu.

**macOS** (Apple Silicon)
- Download `Yleaf.4.1_4.1.5_aarch64.dmg`, open it, and drag **Yleaf 4.1** to your Applications folder.
- Open **Yleaf 4.1** from Applications.

> **"Yleaf 4.1 is damaged and can't be opened" (Gatekeeper warning)**
> macOS blocks apps that are not code-signed with a paid Apple Developer ID certificate.
> To bypass this, run the following command in Terminal **before** opening the DMG:
> ```bash
> xattr -cr ~/Downloads/Yleaf.4.1_4.1.5_aarch64.dmg
> ```
> Then open the DMG and drag the app to Applications. On the first launch you may need to right-click the app → **Open** instead of double-clicking.

**Linux**
- `.AppImage`: download `Yleaf.4.1_4.1.5_amd64.AppImage`, make it executable (`chmod +x`), and run directly.
- `.deb`: install with `sudo dpkg -i Yleaf.4.1_4.1.5_amd64.deb` and launch from your application menu.

No Python, samtools, or other tools need to be installed separately — everything is bundled.

### Portable dashboard (no installation)

If you cannot (or prefer not to) run an installer, download the portable zip for
your platform from the [releases page](https://github.com/genid/Yleaf/releases),
unzip it anywhere (e.g. a USB stick), and run the app directly — no admin rights
and no installation required. Like the installers, everything (Python, samtools,
minimap2, bcftools) is bundled.

**Windows** — `yleaf-dashboard-windows-portable.zip`
- Unzip and double-click `yleaf-dashboard.exe`. Keep `yleaf.exe` next to it (the
  analysis engine the app launches).
- Requires the **Microsoft Edge WebView2 runtime**, which is preinstalled on
  Windows 10/11. On older systems install it once from
  [Microsoft](https://developer.microsoft.com/microsoft-edge/webview2/).
- The app is unsigned, so SmartScreen may show a warning on first launch —
  choose **More info → Run anyway**.

**macOS** (Apple Silicon) — `yleaf-dashboard-macos-portable.zip`
- Unzip and double-click **Yleaf 4.1** to run (no need to move it to Applications).
- The same Gatekeeper note as above applies; if it reports the app is "damaged",
  run `xattr -cr` on the unzipped `.app` first.

**Linux** — `yleaf-dashboard-linux-portable.zip`
- Unzip to get the `.AppImage`, make it executable (`chmod +x`), and run it
  directly.

### UYSD (Y-SNP Database) integration

Each yfull prediction in the result pane shows an embedded haplogroup-frequency
world map sourced from the [Y-SNP Database](https://ysnp.erasmusmc.nl/).
Click the preview to open the full map in your system browser.

You can also contribute your own samples back to UYSD without leaving the
dashboard: the **Submit to UYSD ↗** button opens a panel where you fill in
country / region / publication metadata per sample (with autocompletion against
UYSD's accepted name list and a CSV import/export for bulk editing), then
submit directly.  Requires a UYSD account.

---

## Command-line installation

For users who prefer the command line. Requires Python 3.7+.

### Requirements

    Operating system: Linux, macOS, or Windows (via standalone executable).
    Internet connection: when running for the first time for downloading the reference genome. Alternatively you 
                         can configure your own references.
    Data storage: For installation we recommend a storage capacity of > 8 GB. 

### Option 1: Standalone executable (no Python required)

Download the pre-built binary for your platform from the [releases page](https://github.com/genid/Yleaf/releases):

- `yleaf-linux.tar.gz` — Linux x86-64
- `yleaf-macos.tar.gz` — macOS (Intel and Apple Silicon)
- `yleaf-windows.zip` — Windows x86-64

Extract and run the `yleaf` executable directly. samtools and bcftools are bundled — no external tools needed.

```bash
# Linux / macOS
tar xzf yleaf-linux.tar.gz
./yleaf/yleaf -h

# Windows (PowerShell)
Expand-Archive yleaf-windows.zip
.\yleaf\yleaf.exe -h
```

### Option 2: Conda environment (recommended for source installs)

```bash
# first clone this repository to get the environment_yleaf.yaml
git clone https://github.com/genid/Yleaf.git
cd Yleaf
# create the conda environment — it will be called yleaf
conda env create --file environment_yleaf.yaml
# activate the environment
conda activate yleaf
# pip install the cloned yleaf into your environment. Using the -e flag allows you to modify the config file in your cloned folder
pip install -e .

# verify that Yleaf is installed correctly. You can call this command from any directory on your system
Yleaf -h 
```

### Option 3: Manual install

```bash
# install python 3.7+ and libraries
apt-get install python3
pip3 install pandas numpy
# install external tools
sudo apt-get install minimap2 samtools bcftools
# clone the yleaf repository
git clone https://github.com/genid/Yleaf.git
cd Yleaf
pip install -e .

# verify that Yleaf is installed correctly
Yleaf -h 
```

After installation you can navigate to `yleaf/config.txt` and add custom paths for the reference genomes listed there. This prevents Yleaf from downloading them on first run and allows you to point to an existing reference. Positions are based on hg38, hg19, or T2T (hs1).

## Usage and examples

Here follow some minimal working examples of how to use Yleaf with different input files. There are additional options that can be used to tune how strict Yleaf is as well as options to get private mutations and a graph showing the positioning of predicted haplogroups in the haplogroup tree.

### BAM or CRAM format

    Yleaf -bam file.bam -o bam_output --reference_genome hg38
    Yleaf -cram file.cram -o cram_output --reference_genome hg38

For BAM and CRAM files `-rg` is optional — Yleaf auto-detects the reference build (hg19/hg38/T2T) from the `@SQ` headers. Only specify it explicitly if you want to override the detected build.

### FASTQ (raw reads)

    Yleaf -fastq raw_reads.fastq -o fastq_output --reference_genome hg38

### VCF input

    Yleaf -vcf variants.vcf.gz -o vcf_output --reference_genome hg38
    Yleaf -vcf variants.vcf    -o vcf_output --reference_genome hg38  # plain .vcf also accepted (auto-bgzipped on intake)

Yleaf accepts both bgzipped (`.vcf.gz`) and plain (`.vcf`) input.  Because a
typical VCF only contains variant positions, the state of every other YFull
marker is inferred from the reference genome's chrY sequence so the
prediction quality matches BAM-mode runs.  VCFs from callers that include
`FORMAT/AD` (DeepVariant, GATK HaplotypeCaller, `bcftools mpileup -a AD`)
are preferred; GT-only VCFs are also supported with reduced read-depth
information.

For **targeted-panel** VCFs (e.g. forensic SNP panels), reference inference is
not appropriate: positions outside the panel were never sequenced, so assuming
they match the reference genome can inject spurious calls and even pull the
prediction to the wrong haplogroup.  Use `--no-ref-inference` to predict from the
genotyped markers only:

    Yleaf -vcf panel.vcf.gz -o vcf_output --reference_genome hg38 --no-ref-inference

QC scores are typically lower than a whole-genome run because far fewer markers
are available — this reflects the smaller evidence base, not a problem.

### PLINK / SNP-array

    Yleaf -plink dataset.bed -o plink_output --reference_genome hg38

### With haplogroup tree visualisation and private mutations

    Yleaf -bam file.bam -o bam_output --reference_genome hg38 -dh -p

`-dh` generates a self-contained interactive HTML file with zoomable tree, per-haplogroup tabs, and PDF export.

### Ancient DNA samples

Use the `-aDNA` / `--ancient_DNA` flag when working with ancient DNA. This ignores G>A and C>T mutations, which are common post-mortem deamination artefacts and would otherwise be misinterpreted as derived alleles.

    Yleaf -bam ancient_sample.bam -o output --reference_genome hg38 --ancient_DNA

### Selecting a haplogroup tree

Yleaf supports multiple reference trees. Use the `--tree` flag to select one or more:

| Tree name   | Description                          |
|-------------|--------------------------------------|
| `yfull`     | YFull v14 (default)                  |
| `yfull_v10` | YFull v10.01 (legacy)                |
| `ftdna`     | FTDNA Y-haplotree                    |
| `isogg`     | ISOGG tree                           |

    # Single tree (default is yfull)
    Yleaf -bam file.bam -o output --reference_genome hg38 --tree yfull

    # Multiple trees in one run (single pileup, per-tree prediction)
    Yleaf -bam file.bam -o output --reference_genome hg38 --tree yfull ftdna isogg

### Mixture analysis (forensic)

The `-mix` flag enables forensic mixture deconvolution: Yleaf identifies the contributing Y-haplogroups in a DNA mixture from multiple male donors.

    Yleaf -bam mixture.bam -o output --reference_genome hg38 --tree yfull -mix

Results are written to a `.mix` file per sample. Mixture analysis is tree-aware and supports all reference trees.

### JSON output

Use `--report-json` to write a structured JSON sidecar alongside the normal TSV. The JSON includes the full untruncated marker list, QC scores, and excluded sub-clades for each sample. In multi-tree mode separate files are written per tree (e.g. `report.yfull.json`).

    Yleaf -bam file.bam -o output --reference_genome hg38 --report-json output/report.json

### Using an existing reference genome

By default Yleaf downloads the reference genome on first run. To skip the download, point Yleaf at an existing FASTA with `--ref-fasta`:

    Yleaf -bam file.bam -o output --reference_genome hg38 --ref-fasta /data/hg38.fa

Alternatively, set the `YLEAF_REF_DIR` environment variable to a directory containing files named `hg38.fa` (or `.fasta`/`.fna`) and Yleaf will find the right one automatically. You can also edit `yleaf/config.txt` to set persistent paths for each build.

### Update check

At startup Yleaf asks GitHub whether a newer release exists and, if so, logs a single line
pointing at the releases page. It never modifies your installation — upgrading stays a
deliberate step, which matters because Yleaf is installed in several different ways and
analyses should stay reproducible.

The check is best-effort: it times out after two seconds, caches its answer for 24 hours,
and is silently skipped when the machine is offline. Disable it with `--no-update-check`,
or by setting `YLEAF_NO_UPDATE_CHECK=1` for offline clusters:

    Yleaf -bam file.bam -o output --reference_genome hg38 --no-update-check

## Additional information

For a more comprehensive manual please have a look at the [yleaf_manual](yleaf_manual.pdf).

If you have a bug to report or a question about installation consider sending an email to 
 a.ralf at erasmusmc.nl or create an issue on GitHub.


### References and Supporting Information
A. Ralf, et al., Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data (2018).

https://academic.oup.com/mbe/article/35/5/1291/4922696
