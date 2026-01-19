---
name: bfabricpy
description: Python interface to the B-Fabric laboratory information management system. Use when working with B-Fabric API for CRUD operations on samples, datasets, workunits, resources, or building B-Fabric applications and scripts.
---

# bfabricPy: Python Interface to B-Fabric

Python client library for interacting with the B-Fabric laboratory information management system (LIMS). Provides a general client for CRUD operations, a relational API for simplified read access, and tools for building B-Fabric applications.

**Documentation**: https://fgcz.github.io/bfabricPy/
**GitHub**: https://github.com/fgcz/bfabricPy

## Quick Reference

```python
from bfabric import Bfabric

# Connect to B-Fabric (uses ~/.bfabricpy.yml)
client = Bfabric.connect()

# Read entities
samples = client.read("sample", {"id": [1, 2, 3]})
workunits = client.read("workunit", {"status": "available"})

# Save/update entities
result = client.save("sample", {"name": "New Sample", "containerid": 123})

# Delete entities
client.delete("sample", 456)

# Check existence
exists = client.exists("sample", 123)

# Entity-based API (read-only, with lazy loading)
from bfabric.entities import Sample, Workunit
sample = Sample.find(client, id=123)
workunit = Workunit.find(client, id=456)
```

## Installation

```bash
# For scripts only (isolated environment)
uv tool install bfabric-scripts

# For development (add to project)
# In pyproject.toml:
# dependencies = ["bfabric>=1.16"]

# For app runner
uv tool install bfabric_app_runner
```

## Configuration

Create `~/.bfabricpy.yml`:

```yaml
PRODUCTION:
  login: your_username
  password: your_webservice_password  # NOT your login password
  base_url: https://fgcz-bfabric.uzh.ch/bfabric

# Optional: test environment for integration testing
TEST:
  login: your_username
  password: your_test_password
  base_url: https://fgcz-bfabric-test.uzh.ch/bfabric
```

**Important**: The password is your B-Fabric web service password from your profile page, not your login password.

### Environment Variables

```bash
# Switch between config sections
export BFABRICPY_CONFIG_ENV=TEST

# Override config entirely (JSON format)
export BFABRIC_CONFIG_DATA='{"login": "user", "password": "pass", "base_url": "..."}'
```

## Package Structure

| Package | Purpose |
|---------|---------|
| `bfabric` | Core Python API library |
| `bfabric-scripts` | Command-line utilities |
| `bfabric-app-runner` | Application execution framework |
| `bfabric-asgi-auth` | Authentication middleware |
| `bfabric-rest-proxy` | REST API gateway |

## Client API

### Creating Clients

```python
from bfabric import Bfabric

# Standard connection (respects config hierarchy)
client = Bfabric.connect()

# Specify environment explicitly
client = Bfabric.connect(config_env="TEST")

# Disable config file fallback
client = Bfabric.connect(config_file_env=None)

# For web applications (with token)
client, token_data = Bfabric.connect_webapp(token=request_token)

# Token-based authentication
client = Bfabric.connect_token(token=auth_token)
```

### Configuration Hierarchy

Settings are applied in order (highest to lowest priority):
1. `BFABRIC_CONFIG_DATA` environment variable
2. Config file with explicitly specified environment
3. `BFABRICPY_CONFIG_ENV` environment variable
4. Default environment from `~/.bfabricpy.yml`

### CRUD Operations

```python
from bfabric import Bfabric

client = Bfabric.connect()

# READ - Query entities with filters
samples = client.read(
    endpoint="sample",
    obj={"containerid": 123, "status": "available"},
    max_results=100
)
for sample in samples:
    print(f"Sample {sample['id']}: {sample['name']}")

# READ - Multiple IDs
samples = client.read("sample", {"id": [1, 2, 3, 4, 5]})

# READ - With pagination
all_workunits = client.read(
    "workunit",
    {"applicationid": 456},
    max_results=1000,
    offset=0
)

# SAVE - Create new entity
new_sample = client.save(
    endpoint="sample",
    obj={
        "name": "New Sample",
        "containerid": 123,
        "type": "Biological Sample"
    }
)
sample_id = new_sample[0]["id"]

# SAVE - Update existing entity
client.save(
    endpoint="sample",
    obj={"id": sample_id, "description": "Updated description"},
    method="update"
)

# DELETE - Remove entity
client.delete("sample", sample_id)

# EXISTS - Check if entity exists
if client.exists("sample", sample_id):
    print("Sample exists")

# EXISTS - With additional conditions
if client.exists("sample", sample_id, {"status": "available"}):
    print("Sample exists and is available")
```

### Uploading Resources

```python
from pathlib import Path

# Upload file to workunit (base64 encoded)
client.upload_resource(
    resource_name="results.txt",
    content=Path("results.txt").read_bytes(),
    workunit_id=12345
)
```

### Context Manager for Auth Switching

```python
# Temporarily switch authentication
with client.with_auth(login="other_user", password="other_pass"):
    # Operations run as other_user
    samples = client.read("sample", {"containerid": 123})
# Back to original auth
```

## Entity API (Read-Only)

The `bfabric.entities` module provides a read-only API with lazy loading for associated entities.

### Basic Usage

```python
from bfabric import Bfabric
from bfabric.entities import Sample, Workunit, Resource, Dataset

client = Bfabric.connect()

# Find single entity
sample = Sample.find(client, id=123)
print(f"Sample: {sample.name}")

# Find multiple entities (returns dict[id, entity])
samples = Sample.find_all(client, ids=[1, 2, 3])
for sample_id, sample in samples.items():
    print(f"{sample_id}: {sample.name}")

# Find by filter
available_workunits = Workunit.find_by(
    client,
    {"status": "available", "applicationid": 456}
)
```

### Entity Relationships (Lazy Loading)

```python
# HasOne relationship - loads on access
workunit = Workunit.find(client, id=123)
application = workunit.application  # Lazy loaded
print(f"Application: {application.name}")

# HasMany relationship - loads on access
container = Container.find(client, id=456)
samples = container.samples  # Lazy loaded list
for sample in samples:
    print(f"Sample: {sample.name}")

# Optional relationships return None or empty list
resource = Resource.find(client, id=789)
dataset = resource.dataset  # None if not linked
```

### Available Entities

Common entity types (each with appropriate ENDPOINT):
- `Sample` - Biological samples
- `Workunit` - Processing units
- `Resource` - Files and outputs
- `Dataset` - Data collections
- `Container` - Sample containers (plates, tubes)
- `Application` - B-Fabric applications
- `Project` - Research projects
- `User` - System users

## Command-Line Scripts

After installing `bfabric-scripts`:

```bash
# List available commands
bfabric --help

# Read operations
bfabric read sample --id 123
bfabric read workunit --filter '{"status": "available"}'

# The scripts log active user and URL for verification
```

## App Runner

The `bfabric-app-runner` enables building B-Fabric applications with standardized input/output handling.

### Installation

```bash
uv tool install bfabric_app_runner
```

### WorkunitDefinition

```python
from bfabric_app_runner import WorkunitDefinition

# Create workunit definition
workunit_def = WorkunitDefinition(
    workunit_id=12345,
    execution_params={"param1": "value1"},
    registration_params={"output_type": "dataset"}
)

# Convert to/from YAML
yaml_str = workunit_def.to_yaml()
workunit_def = WorkunitDefinition.from_yaml(yaml_str)
```

### Input Specification

Define how input files are validated and prepared:

```yaml
inputs:
  - name: raw_data
    type: resource
    required: true
    validation:
      extension: [.raw, .mzML]
```

### Output Specification

Define how results are registered:

```yaml
outputs:
  - name: results
    type: resource
    copy_to: workunit
  - name: summary
    type: dataset
    save: true
```

### App Specification

Define application execution:

```yaml
app:
  name: my_analysis
  version: "1.0.0"

  # Direct command execution
  command: python analyze.py --input {input_file} --output {output_dir}

  # Or Docker execution
  docker:
    image: myregistry/myapp:1.0.0
    command: /app/run.sh
```

### Template Variables

Use variable substitution in commands:

```yaml
command: >
  python process.py
    --input {inputs.raw_data.path}
    --output {outputs.results.path}
    --workunit {workunit_id}
    --param {execution_params.threshold}
```

## Common Patterns

### Pattern 1: Process All Samples in Container

```python
from bfabric import Bfabric
from bfabric.entities import Container, Sample

client = Bfabric.connect()

container = Container.find(client, id=12345)
for sample in container.samples:
    print(f"Processing {sample.name}")
    # Process sample...

    # Update sample status
    client.save("sample", {
        "id": sample.id,
        "status": "processed"
    }, method="update")
```

### Pattern 2: Create Workunit with Resources

```python
from bfabric import Bfabric
from pathlib import Path

client = Bfabric.connect()

# Create workunit
workunit = client.save("workunit", {
    "name": "Analysis Results",
    "applicationid": 123,
    "containerid": 456,
    "status": "available"
})[0]

workunit_id = workunit["id"]

# Upload result files
for result_file in Path("results/").glob("*.txt"):
    client.upload_resource(
        resource_name=result_file.name,
        content=result_file.read_bytes(),
        workunit_id=workunit_id
    )

print(f"Created workunit {workunit_id} with results")
```

### Pattern 3: Query and Export Data

```python
from bfabric import Bfabric
import pandas as pd

client = Bfabric.connect()

# Query samples with specific criteria
samples = client.read("sample", {
    "containerid": 123,
    "type": "Biological Sample",
    "status": "available"
}, max_results=1000)

# Convert to DataFrame
df = pd.DataFrame(samples)
df.to_csv("samples_export.csv", index=False)
```

### Pattern 4: Batch Processing with Error Handling

```python
from bfabric import Bfabric
from loguru import logger

client = Bfabric.connect()

sample_ids = [1, 2, 3, 4, 5]
results = []

for sample_id in sample_ids:
    try:
        sample = client.read("sample", {"id": sample_id})
        if sample:
            # Process sample
            result = process_sample(sample[0])
            results.append(result)
            logger.info("Processed sample {}", sample_id)
        else:
            logger.warning("Sample {} not found", sample_id)
    except Exception as e:
        logger.error("Error processing sample {}: {}", sample_id, e)
        continue

logger.info("Processed {} of {} samples", len(results), len(sample_ids))
```

### Pattern 5: Web Application Integration

```python
from bfabric import Bfabric

def handle_webapp_request(token: str):
    """Handle request from B-Fabric web application."""
    # Connect with webapp token
    client, token_data = Bfabric.connect_webapp(token=token)

    # Access entity info from token
    entity_class = token_data.entity_class
    entity_id = token_data.entity_id

    # Read the relevant entity
    entity = client.read(entity_class, {"id": entity_id})

    return process_entity(entity[0])
```

## Best Practices

### Configuration
- Never commit credentials to version control
- Use `~/.bfabricpy.yml` for personal credentials
- Use environment variables for CI/CD and production
- Test against TEST environment before production

### Error Handling
- Always check if read operations return results
- Use try/except for network operations
- Log operations for debugging

### Performance
- Use `max_results` to limit large queries
- Batch operations when possible
- Use entity API lazy loading to avoid unnecessary requests

### Development Workflow
- Install `bfabric-scripts` for quick CLI access
- Use TEST environment for development
- Verify user/URL from script logs before operations

## Resources

- **Documentation**: https://fgcz.github.io/bfabricPy/
- **App Runner Docs**: https://fgcz.github.io/bfabricPy/app_runner/
- **GitHub Repository**: https://github.com/fgcz/bfabricPy
- **Issues**: https://github.com/fgcz/bfabricPy/issues

### Citation

When using bfabricPy in research, please cite:

> Bridging data management platforms and visualization tools to enable ad-hoc and smart analytics in life sciences.
> Journal of Integrative Bioinformatics, 2022. doi: 10.1515/jib-2022-0031

## Common Issues

### Authentication Errors

```python
# Verify your credentials
from bfabric import Bfabric

try:
    client = Bfabric.connect()
    # Test with a simple read
    result = client.read("user", {"login": "your_username"})
    print(f"Connected as: {result[0]['login']}")
except Exception as e:
    print(f"Authentication failed: {e}")
    # Check: Is password the WEB SERVICE password (not login)?
    # Check: Is base_url correct?
```

### Entity Not Found

```python
# Always check for empty results
samples = client.read("sample", {"id": 99999})
if not samples:
    print("Sample not found")
else:
    sample = samples[0]
```

### Permission Denied

```python
# Check if you have access to the container/project
# Some operations require specific roles
try:
    client.save("sample", {...})
except Exception as e:
    if "permission" in str(e).lower():
        print("You don't have write access to this container")
```

### Connection Timeout

```python
# For long-running operations, consider pagination
all_results = []
offset = 0
batch_size = 100

while True:
    batch = client.read(
        "sample",
        {"containerid": 123},
        max_results=batch_size,
        offset=offset
    )
    if not batch:
        break
    all_results.extend(batch)
    offset += batch_size
```
