# Convert `raw` files to mzml

```bash
# Create lists of .raw mass spec files split into groups of 8
cd ProtPipe
find /data/CARD_ARDIS/Neurite_massspec_raw/ -name '*.raw'  | \
split -l 8 -d --additional-suffix '.txt' - rawfiles_

sbatch src/pwiz-convert.sh --list rawfiles_00.txt
sbatch src/pwiz-convert.sh --list rawfiles_01.txt
sbatch src/pwiz-convert.sh --list rawfiles_02.txt
sbatch src/pwiz-convert.sh --list rawfiles_03.txt
sbatch src/pwiz-convert.sh --list rawfiles_04.txt
sbatch src/pwiz-convert.sh --list rawfiles_05.txt
sbatch src/pwiz-convert.sh --list rawfiles_06.txt
sbatch src/pwiz-convert.sh --list rawfiles_07.txt
sbatch src/pwiz-convert.sh --list rawfiles_08.txt
```

# Creat config file
`config.txt` points to the `mzML` directory with all mass spec files to be analyzed

# Run DIA-NN
```bash
sbatch src/diann.sh --cfg config.txt
```

Rename design matrix with `sample-rename.R`


# ProtPipe analysis of knockdown experiment
```bash
module load R/4.3
Rscript src/basic_analysis_KD.R \
    --pgfile KD.pg_matrix.tsv \
    --design design_matrix_KD.csv \
    --out Proteomics_KD
```

# ProtPipe analysis of WT Neurite vs Soma
```bash
Rscript src/basic_analysis_WT.R \
    --pgfile WT.pg_matrix.tsv \
    --design design_matrix_02.csv \
    --out Proteomics_WT
```