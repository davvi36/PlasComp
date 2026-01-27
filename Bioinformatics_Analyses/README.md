Supplementary data (`data/gkab859_supplementary_data.xlsx` and `data/gkab859_supplementary_spacers.tsv`) were retrieved from published article: *"CRISPR-Cas systems are widespread accessory elements across bacterial and archaeal plasmids"* (https://doi.org/10.1093/nar/gkab859).
Due to its size, `data/gkab859_supplementary_spacers.tsv` is provided as a .gz and should be extracted before use.

TA nucleotide (`data/tadb.fna`) and protein (`data/tadb.faa`) sequences were retrieved from the TADB3.0 database (https://bioinfo-mml.sjtu.edu.cn/TADB3), and the metadata (`data/tadb_meta.tsv`) collated by hand from the database.

Install pixi if you don't have it yet:

```sh
# curl -fsSL https://pixi.sh/install.sh | sh
```

Install other dependencies:

```sh
pixi install
```

Get plasmid ids:

```sh
pixi run duckdb -c "COPY (SELECT Acc FROM read_xlsx('data/gkab859_supplementary_data.xlsx', sheet='all_plasmids', header=true) WHERE Derep = 'TRUE' AND Acc IS NOT NULL) TO 'plasmid_ids.txt' (HEADER 0);"
```

Download plasmid sequences from NCBI:

```sh
pixi run 'ids=$(paste -sd, data/plasmid_ids.txt); efetch -db nucleotide -id "$ids" -format fasta > data/plasmid_seq.fna'
```

Separate plasmid sequences into individual files:

```sh
mkdir -p data/plasmids_separated/
pixi run seqkit split --by-id data/plasmid_seq.fna --out-dir data/plasmids_separated/ --by-id-prefix '' --threads 10
```

Predict plasmid ORFs:

```sh
mkdir -p data/plasmid_orfs/
pixi run parallel -j 14 --progress "prodigal -i data/plasmids_separated/{}.fna -a data/plasmid_orfs/{}.faa -p meta -f gff -o data/plasmid_orfs/{}.gff" :::: data/plasmid_ids.txt
```

Concatenate plasmid protein sequences and GFF:

```sh
cat data/plasmid_orfs/*.faa > data/plasmid_seq.faa
cat data/plasmid_orfs/*.gff > data/plasmid_seq.gff
```

Make BLAST databases for plasmid sequences:

```sh
mkdir -p data/blastdb
pixi run makeblastdb -in "data/plasmid_seq.faa" -dbtype prot
mv data/plasmid_seq.faa.* data/blastdb/
pixi run makeblastdb -in "data/plasmid_seq.fna" -dbtype nucl
mv data/plasmid_seq.fna.* data/blastdb/
```

Search for TA systems in plasmid ORFs:

```sh
mkdir -p "data/ta_blast/"
pixi run blastp -query "data/tadb.faa" -db "data/blastdb/plasmid_seq.faa" -evalue 0.05 -outfmt 6 -num_threads 10 -out "data/ta_blast/plasmid_faa.m8"
pixi run blastn -query "data/tadb.fna" -db "data/blastdb/plasmid_seq.fna" -evalue 0.05 -outfmt 6 -num_threads 10 -out "data/ta_blast/plasmid_fna.m8"
```

Run `analysis.R` interactively to generate tables and plots.
