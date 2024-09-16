# 
```bash
# Clone ProtPipe repository from GitHub
git clone https://github.com/NIH-CARD/ProtPipe.git


# Retrieve UniProt human proteome
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz

# Unzip into ProtPipe
gunzip -c UP000005640_9606.fasta.gz > ProtPipe/UP000005640_9606.fasta && \
rm UP000005640_9606.fasta.gz
```

See `Proteomics` directory for more steps