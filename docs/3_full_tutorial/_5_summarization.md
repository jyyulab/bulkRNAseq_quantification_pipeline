---
title: Summarization
layout: default
nav_order: 5
parent: Full Tutorial
---


## Quantification Summary

For each single sample, we will generate a QC report from the quantification results by Salmon and RSEM. There are five QC metrics in each report:

* Alignment statistics: showing read counts and mapping rates et. al.;
* Quantification statistics: showing the number of genes/transcripts identified by two methods and the correlations of quantification results of them;
* Biotype distribution: showiing the composition of types of identified transcripts and genes;
* Quantification accuracy: correlations of gene expression measurements by Salmon and RSEM;
* Genebody coverage statistics: showing if the RNA samples were degraded or not.

![image-20230901163554962](/Users/qpan/Library/Application Support/typora-user-images/image-20230901163554962.png)

```sh
Rscript -e "rmarkdown::render(input = '/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/bin/summaryIndividual.Rmd', clean = TRUE, quiet = F, output_format = 'html_document', output_file = 'quantSummary.html', output_dir = '/your_path/quantSummary', params = list(sampleName = 'sample1', quant_dir = '/your_path'))"
```

As shown above, we built a R markdown script, **summaryIndividual.Rmd**, to generate this report from the outputs of quantification pipelines. The only argument you need to pay attention to is the **quant_dir**. The script will call these files to generate the QC report:

* quant_dir/preProcessing/adapterTrimming.json
* quant_dir/quantRSEM/quant.isoforms.results
* quant_dir/quantRSEM/quant.genes.results
* quant_dir/quantRSEM/quant.stat/quant.cnt
* quant_dir/quantRSEM/geneCoverage.txt
* quant_dir/quantSalmon/quant.sf
* quant_dir/quantSalmon/quant.genes.sf

As for the outputs, the script will generate a .html file named quantSummary.html which summarizes all five QC metrics. And a few .txt files and .pdf files will be generated as well to provide the source data of the html QC report. These files can be also used to generate the multi-sample QC report if multiple samples were sequenced.