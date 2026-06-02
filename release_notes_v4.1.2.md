## Yleaf v4.1.2

- ISOGG `~` (approximate-placement) markers no longer count toward any QC score; regenerated `isogg_positions_{hg19,hg38,t2t}.txt` preserves the `~` from the source so e.g. R1b1a1b1a1a1~ markers (FGC36477/8/9) no longer drag R1b1a1b1a1a1's QC2 below threshold.
- Dashboard: Browse file accepts multiple inputs in one shot; they are staged into a temp dir (with their `.bai`/`.crai`/`.csi` indexes) and processed by the directory-mode pipeline.
- Dashboard: directories containing multiple eligible input types (e.g. BAM + VCF) now fan out one job per type into `<output_dir>/<type>/`, instead of silently picking only the highest-priority type.
