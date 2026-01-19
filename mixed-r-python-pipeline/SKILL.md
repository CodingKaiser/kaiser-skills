---
name: mixed-r-python-pipeline
description: Expert guidance for turning existing Snakemake analysis projects (with Snakefile + src/) into pip-installable CLI tools that can deploy the workflow to new projects. Use when the user has a working analysis project with Snakemake and R/Python scripts and wants to make it reusable/deployable. Creates CLI with init, validate, update, clean, info commands. Reference implementation at /Users/wolski/projects/ptm-pipeline.
---

# Turn Existing Snakemake Projects into Deployable CLI Tools

Take an existing analysis project with a working Snakefile and src/ directory and wrap it into a pip-installable CLI tool.

## Starting Point

You already have a working project:
```
my-analysis-project/
├── Snakefile              ← Already exists and works
├── helpers.py             ← Already exists
├── src/                   ← Already exists
│   ├── analysis.Rmd
│   ├── process.R
│   └── integrate.py
├── config.yaml            ← Project-specific config
└── data/                  ← Project data
```

## Goal

Create a CLI tool that can deploy this workflow to new projects:
```bash
my-pipeline init /path/to/new/project
# → Discovers data, generates config, copies Snakefile + src/
```

## Reference Implementation

**ptm-pipeline**: `/Users/wolski/projects/ptm-pipeline`

Key files to study:
- `src/ptm_pipeline/cli.py` - Click commands
- `src/ptm_pipeline/discover.py` - Glob-based data discovery
- `src/ptm_pipeline/init.py` - Template copying logic
- `src/ptm_pipeline/config.py` - YAML generation
- `template/` - The original Snakefile + src/ moved here

## Transformation Steps

### Step 1: Create CLI Project Structure

```bash
mkdir my-pipeline && cd my-pipeline
uv init

mkdir -p src/my_pipeline template
```

### Step 2: Move Existing Workflow to template/

```bash
# From your existing working project:
cp Snakefile ../my-pipeline/template/
cp helpers.py ../my-pipeline/template/
cp -r src/ ../my-pipeline/template/src/
cp Makefile ../my-pipeline/template/  # if exists
```

### Step 3: Create CLI Modules

Create in `src/my_pipeline/`:

| Module | Purpose |
|--------|---------|
| `cli.py` | Click commands: init, validate, update, clean, info |
| `discover.py` | Auto-detect data folders with glob patterns |
| `init.py` | Copy template, generate config |
| `config.py` | YAML config generation |
| `validate.py` | Check files, commands, R packages |
| `clean.py` | Remove pipeline files |

### Step 4: Configure pyproject.toml

```toml
[project]
name = "my-pipeline"
version = "0.1.0"
dependencies = ["click>=8.0", "pyyaml>=6.0", "rich>=13.0"]

[project.scripts]
my-pipeline = "my_pipeline.cli:main"

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

# CRITICAL: This makes template/ installable
[tool.hatch.build.targets.wheel.shared-data]
"template" = "share/my-pipeline/template"
```

## Key Patterns

### Template Path Resolution

```python
def get_template_path() -> Path:
    # 1. Development mode
    dev_path = Path(__file__).parent.parent.parent / "template"
    if dev_path.exists():
        return dev_path

    # 2. Installed mode
    installed_path = Path(sys.prefix) / "share" / "my-pipeline" / "template"
    if installed_path.exists():
        return installed_path

    raise FileNotFoundError("Could not locate template directory")
```

### Data Discovery

```python
def find_folders(directory: Path, patterns: List[str]) -> Dict[str, List[Path]]:
    results = {}
    for pattern in patterns:
        matches = sorted(directory.glob(pattern))
        if matches:
            key = pattern.replace("*", "").strip("_")
            results[key] = matches
    return results

# Define patterns for YOUR data types:
DATA_PATTERNS = ["DEA_*", "results_*", "input_*"]
```

### Copy Template to Target

```python
def copy_template_files(target_dir: Path, template_path: Path):
    for item in ["Snakefile", "helpers.py", "Makefile"]:
        src = template_path / item
        if src.exists():
            shutil.copy2(src, target_dir / item)

    src_dir = template_path / "src"
    if src_dir.exists():
        target_src = target_dir / "src"
        if target_src.exists():
            shutil.rmtree(target_src)
        shutil.copytree(src_dir, target_src)
```

## Result

```
my-pipeline/                          # The CLI tool
├── pyproject.toml
├── src/my_pipeline/
│   ├── cli.py
│   ├── discover.py
│   ├── init.py
│   ├── config.py
│   ├── validate.py
│   └── clean.py
└── template/                         # Your original workflow
    ├── Snakefile
    ├── helpers.py
    └── src/
        ├── analysis.Rmd
        └── ...

# Usage:
uv run my-pipeline init /path/to/new/project
uv run my-pipeline validate /path/to/new/project
cd /path/to/new/project && snakemake -s Snakefile -c1 all
```

## CLI Commands

| Command | Purpose |
|---------|---------|
| `init [DIR]` | Discover data, prompt settings, generate config, copy template |
| `validate [DIR]` | Check files, commands, R packages exist |
| `update [DIR]` | Re-copy template files, preserve config |
| `clean [DIR]` | Remove pipeline files, keep user data |
| `info [DIR]` | Show discovered data folders (debugging) |
