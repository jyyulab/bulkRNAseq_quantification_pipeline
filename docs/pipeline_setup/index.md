---
title: Pipeline Setup
layout: default
nav_order: 2
permalink: /docs/pipeline_setup/index
---

### Before you start

There are two parts to set up this pipeline:

- **Software installation**: to install all tools required by this pipeline.
- **Reference preparation**: to generate reference files used in this pipeline.

For Yu Lab members, we have set this pipeline up on HPC:

- All the tools have been compiled in one conda environment, which can be launched by:

  ``` bash
  module load conda3/202402
  conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2025
  ```

- We have generated the files for four reference genomes: hg38, hg19, mm39 and mm10.


  | Genome            | GENCODE release | Release date | Ensembl release | Path                                                         |
  | ----------------- | --------------- | ------------ | --------------- | ------------------------------------------------------------ |
  | hg38 (GRCh38.p14) | v48             | 05.2025      | v114            | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release48 |
  | hg19 (GRCh37.p13) | v48lift37*      | 05.2025      | v114            | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release48 |
  | mm39 (GRCm39)     | vM37            | 05.2025      | v114            | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM37 |
  | mm10 (GRCm38.p6)  | vM25            | 04.2020**    | v100            | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25 |

  *: The updates for the hg19/GRCh37 genome assembly have stopped in 2013. However, gene annotation continue to be updated by mapping the comprehensive gene annotations originally created for the GRCh38/hg38 reference chromosomes onto GRCh37 primary assembly using [gencode-backmap](https://github.com/diekhans/gencode-backmap) .

  **: The updates for the mm10/GRCm38 genome assembly and gene annotation have stopped in 2019.

You should run this tutorial only when you want to set up this pipeline locally or complie other reference genome / gene annotation.