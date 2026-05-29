## Yleaf v4.1.0

- Dashboard now integrates with the [Y-SNP Database (UYSD)](https://ysnp.erasmusmc.nl/): embedded haplogroup-frequency maps in the result pane, plus a built-in submission flow.
- VCF input now produces real predictions on whole-genome inputs: reference-based ancestral inference restored, plain `.vcf` files accepted (resolves issue #54).
- `--ref-fasta` is validated against `--reference_genome` and raises a clear error on build mismatch.
- Ancient-DNA position files regenerated as main minus `C→T`/`G→A` rows; marker counts roughly double across all builds.
- T2T marker names and chrY reference restored — both were shipping as empty/blank placeholders.
- Python 3.7 compatibility restored for the bundled sidecar.
