---
name: snakemake-compact
description: Compact Snakemake workflow patterns. Keep rules short, complex logic in Python modules.
---

# Snakemake Essentials

## File Naming
- Snakefiles: `workflow.smk`, `rules.smk`, `Snakefile.smk`
- Helper modules: `helpers.py` next to Snakefile
- Config: `config.yaml`

## Basic Rule (quote paths; named IO)
Use `:q` so file paths with spaces don't break.

Always show details
```python
rule align:
    input:
        reads="data/{sample}.fastq",
        index="ref/genome.idx"
    output:
        bam="results/{sample}.bam"
    log:
        "logs/align/{sample}.log"
    benchmark:
        "bench/align/{sample}.tsv"
    shell:
        "aligner -i {input.reads:q} -x {input.index:q} -o {output.bam:q} 2> {log:q}"
```

## Complex Logic in Python Module
```python
# workflow.smk
from helpers import run_analysis

rule analyze:
    input:
        data="data/{sample}.csv"
    output:
        report="results/{sample}_analysis.json"
    run:
        run_analysis(input.data, output.report)
```

```python
# helpers.py (next to workflow.smk)
import json
from pathlib import Path

def run_analysis(input_file: str, output_file: str) -> None:
    """Complex analysis logic lives here."""
    data = Path(input_file).read_text()
    result = {"status": "ok", "lines": len(data.splitlines())}
    Path(output_file).write_text(json.dumps(result))
```

## Core Patterns (keep Snakefile compact)
```python
# Keep Snakefile minimal; include rule files
include: "rules/qc.smk"
include: "rules/align.smk"

# Override ambiguities and force local steps
ruleorder: fast_align > align
localrules: prep_refs

# Tools that emit directories
rule assemble:
    output:
        outdir=directory("results/{sample}/assembly")

# Sentinel when outputs are many/variable
rule done:
    input:
        outdir="results/{sample}/assembly"
    output:
        done=touch("results/{sample}/assembly.done")
```

## Wildcards & Expansion
```python
SAMPLES = ["A", "B", "C"]

rule all:
    input:
        expand("results/{sample}.bam", sample=SAMPLES)

rule process:
    input:
        reads="data/{sample}.txt"
    output:
        bam="results/{sample}.bam"
    shell:
        "process {input.reads:q} > {output.bam:q}"
```

## Config
```yaml
# config.yaml
samples: ["A", "B", "C"]
threads: 8
reference: "ref/genome.fa"
```

```python
configfile: "config.yaml"

rule align:
    threads: config["threads"]
    input:
        reads="data/{sample}.fastq",
        ref=config["reference"]
    output:
        bam="results/{sample}.bam"
    shell:
        "aligner -i {input.reads:q} -x {input.ref:q} -o {output.bam:q}"
```

## Params
```python
rule filter:
    input:
        vcf="data/{sample}.vcf"
    output:
        vcf="filtered/{sample}.vcf"
    params:
        qual=30,
        extra=lambda wc: f"--sample {wc.sample}"
    shell:
        "bcftools filter -q {params.qual} {params.extra} {input.vcf:q} > {output.vcf:q}"
```

## Temp & Protected Files
```python
rule step1:
    output:
        tmp=temp("intermediate/{sample}.tmp")  # Auto-deleted

rule final:
    input:
        tmp="intermediate/{sample}.tmp"
    output:
        final=protected("results/{sample}.final")  # Read-only
```

## CLI
```bash
snakemake -n -j 1               # Dry-run
snakemake -j 4                  # Run with 4 cores
snakemake --dag -j 1 | dot -Tpng > dag.png  # Visualize DAG
snakemake --lint -j 1           # Check workflow
snakemake -F -j 1               # Force re-run all
snakemake target.txt -j 1       # Build specific target
```

## Key Principles
1. **Short rules**: Import functions from `helpers.py`, call in `run:` block. (Replaces `script:` directive for unified `uv`/pip envs).
2. **Use `.smk` extension** for Snakefiles
3. **One input/output per line** for readability
4. **Named inputs/outputs** over positional: `input: reads="..."` not `input: "..."`
5. **Config over hardcoding**: Put paths/params in `config.yaml`
6. **KISS**: Simple rules, modular helper functions
7. **Deterministic outputs**: Avoid timestamps in output filenames
