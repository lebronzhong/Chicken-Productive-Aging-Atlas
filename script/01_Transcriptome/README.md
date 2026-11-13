# Transcriptome reconstruction and annotation

The 5 principal analysis entries are listed below. Grouped entries link to their component scripts; other files provide supporting functions.

| Analysis | Main code |
|---|---|
| Read processing | [process_reads.sh](01_Read_processing/process_reads.sh) |
| Transcript reconstruction and annotation merge | [reconstruct_transcripts.sh](02_Assembly/reconstruct_transcripts.sh); [merge_annotations.sh](02_Assembly/merge_annotations.sh) |
| ORF annotation | [annotate_orfs.sh](03_ORF_annotation/annotate_orfs.sh) |
| Six-resource comparison | [query_oriented_six_resources.py](04_Reference_comparison/query_oriented_six_resources.py) |
| Short-read validation | [validation.py](05_Short_read_validation/validation.py) |

## Parameterized invocation

Replace each `{...}` placeholder with your input value.

```bash
bash script/01_Transcriptome/01_Read_processing/process_reads.sh "{sample_directory}" "{genome_fasta}"
bash script/01_Transcriptome/02_Assembly/reconstruct_transcripts.sh "{bam}" "{genome_fasta}" "{reference_gtf}" "{output_dir}"
bash script/01_Transcriptome/02_Assembly/merge_annotations.sh "{reference_gtf}" "{gtf_list}" "{output_prefix}"
```

For programs that locate files relative to their script, prepare the working directory with `tools/prepare_workspace.py` and supply the [required inputs](../../docs/INPUTS.md). Install the external programs in the environment for the selected analysis.

`annotate_orfs.sh` consolidates the original ORF steps. Set `TRANSDECODER_BIN` and choose `sequence`, `long_orfs`, `homology`, `domains`, `predict` or `project`. The script checks the required arguments for each stage. Reference comparison configurations are supplied beside the projection programs.
