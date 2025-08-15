# dynaest

## description

- `dynaest` is a tool for estimating the metabolic labeling scRNA-seq data.
- `dynaest` inherits core algorithms from [dynast-release](https://github.com/aristoteleo/dynast-release.git), and maintains compatibility with `celescope` output formats.


## installation

```
git clone https://github.com/pigraul/dynaest.git
cd dynaest
pip install .
```


## usage

- step1: `dedup`
  - `dedup` is used to remove UMI duplicate reads in a bam file.
  - set `--input` to the bam file from conversion step.
  - set `--output` to the output file.
  - For large datasets, set `-t 1` to reduce memory usage.

  `dynaest dedup --input ../bams/test.PosTag.bam --output count `


- step2: `estimate`
  - `estimate` is used to perform statistical estimation through expectation maximization (EM) and Bayesian inference.
  - set `--gtf` to the gtf file. 
  - set `-o` to the output file.
  - set the output dir from `dedup` step

  ```
  dynaest estimate  \
   -o estimate \
   --gtf /Database/genome/mus_musculus/Mus_musculus.GRCm38.92.chr.gtf \
   ../count
  ```
