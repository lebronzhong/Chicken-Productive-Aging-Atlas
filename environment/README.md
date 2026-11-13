# Software and environments

Use separate environments for ONT processing, expression analysis, molecular QTLs, allelic analysis, sequence modeling and LDSC. Install the dependencies needed for the selected module.

| Stage | Software |
|---|---|
| ONT reads | Guppy 5.0.16; NanoFilt 2.8.0; SeqKit 0.12.0; Pychopper 2.4.0; minimap2 2.26-r1175; SAMtools 1.2 |
| Transcript reconstruction / requantification | IsoQuant 3.3.1 / 3.5.1; gffcompare 0.12.9 |
| Reference projection | Liftoff 1.6.3; minimap2 2.28-r1209 |
| Short-read validation | STAR 2.7.11b; featureCounts 2.0.8 |
| ORF annotation | TransDecoder 5.7.1; DIAMOND 2.1.8.162; HMMER 3.4; SQANTI3 |
| Specificity and functional annotation | tspex 0.6.2; org.Gg.eg.db 3.22.0; GO.db 3.22.0; MEME/AME 5.5.7 |
| Expression and transcript usage | DESeq2; edgeR; SUPPA |
| WGS and allelic analysis | fastp 0.23.4; BWA 0.7.17; GATK 4.1.4.1 and 4.3.0.0; BCFtools 1.21; LORALS; minimap2 2.31 for transcript alignment |
| Molecular QTLs | Salmon; tensorQTL 1.0.9/1.0.10; coloc 5.2.3; SMR 1.03 |
| Chicken GWAS set tests | R 4.5.1; BEDTools 2.31.1 |
| Human trait enrichment | LDSC 1.0.1 and a compatible Python environment |
| DNA sequence model | Python 3.10.6; TensorFlow 2.15.1; CUDA 12.2; Baskerville TF checkout `b06950be1be704ab09ab23ffb6f5102d9e2a63a0` |

Python and R package lists are provided in `python-analysis.txt`, `python-model.txt` and `R-packages.txt`. LORALS reference preparation requires pysam; SMR input preparation requires pandas and pyarrow. Shell scripts use Bash and GNU command-line tools.

Activate the appropriate environment before running a module. Configure executable paths, scheduler resources and input manifests for the computing system. The DNA model settings are in `script/05_DNA_sequence_model/params_medium_33.json`.
