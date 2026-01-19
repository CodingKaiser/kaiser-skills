# Python Anti-Patterns and Fixes

Common Python mistakes and their corrections.

## 1. Mutable Default Arguments

**Anti-pattern:**
```python
def add_item(item, items=[]):  # WRONG
    items.append(item)
    return items
```

**Why it's wrong:** The list is created once when the function is defined, not each time it's called.

**Fix:**
```python
def add_item(item, items=None):
    if items is None:
        items = []
    items.append(item)
    return items
```

## 2. Bare Except Clauses

**Anti-pattern:**
```python
try:
    risky_operation()
except:  # WRONG - catches everything, including KeyboardInterrupt
    handle_error()
```

**Fix:**
```python
try:
    risky_operation()
except Exception as e:  # Or specific exception types
    logger.error(f"Operation failed: {e}")
    handle_error()
```

## 3. Using == for None Comparisons

**Anti-pattern:**
```python
if value == None:  # WRONG
    ...
```

**Fix:**
```python
if value is None:
    ...
```

**Why:** `is` checks identity, `==` checks equality. `None` is a singleton.

## 4. Comparing Boolean Values Explicitly

**Anti-pattern:**
```python
if flag == True:  # WRONG
    ...
if len(items) > 0:  # WRONG
    ...
```

**Fix:**
```python
if flag:
    ...
if items:
    ...
```

## 5. Not Using Context Managers for Files

**Anti-pattern:**
```python
f = open("file.txt")  # WRONG - file may not close if error occurs
data = f.read()
f.close()
```

**Fix:**
```python
with open("file.txt") as f:
    data = f.read()
```

## 6. String Concatenation in Loops

**Anti-pattern:**
```python
result = ""
for item in items:
    result += str(item)  # WRONG - creates new string each iteration
```

**Fix:**
```python
result = "".join(str(item) for item in items)
```

## 7. Modifying List While Iterating

**Anti-pattern:**
```python
for item in items:
    if should_remove(item):
        items.remove(item)  # WRONG - skips elements
```

**Fix:**
```python
items = [item for item in items if not should_remove(item)]
# Or
items[:] = [item for item in items if not should_remove(item)]
```

## 8. Using eval() or exec()

**Anti-pattern:**
```python
user_input = get_user_input()
result = eval(user_input)  # WRONG - major security risk
```

**Fix:**
```python
import ast
result = ast.literal_eval(user_input)  # Only evaluates literals
```

## 9. Not Using enumerate()

**Anti-pattern:**
```python
i = 0
for item in items:
    print(f"{i}: {item}")
    i += 1
```

**Fix:**
```python
for i, item in enumerate(items):
    print(f"{i}: {item}")
```

## 10. Creating Empty Lists/Dicts Unnecessarily

**Anti-pattern:**
```python
items = []
items.append(1)
items.append(2)
items.append(3)
```

**Fix:**
```python
items = [1, 2, 3]
```

## 11. Not Using dict.get() with Defaults

**Anti-pattern:**
```python
if key in my_dict:
    value = my_dict[key]
else:
    value = default
```

**Fix:**
```python
value = my_dict.get(key, default)
```

## 12. Using range(len()) Instead of enumerate()

**Anti-pattern:**
```python
for i in range(len(items)):
    item = items[i]
    print(f"{i}: {item}")
```

**Fix:**
```python
for i, item in enumerate(items):
    print(f"{i}: {item}")
```

## 13. Not Using Collections Module

**Anti-pattern:**
```python
word_counts = {}
for word in words:
    if word in word_counts:
        word_counts[word] += 1
    else:
        word_counts[word] = 1
```

**Fix:**
```python
from collections import Counter
word_counts = Counter(words)
```

## 14. Not Using defaultdict

**Anti-pattern:**
```python
groups = {}
for item in items:
    key = get_key(item)
    if key not in groups:
        groups[key] = []
    groups[key].append(item)
```

**Fix:**
```python
from collections import defaultdict
groups = defaultdict(list)
for item in items:
    key = get_key(item)
    groups[key].append(item)
```

## 15. Overly Complex Comprehensions

**Anti-pattern:**
```python
result = [
    transform(x)
    for x in items
    if condition1(x)
    if condition2(x)
    if condition3(x)
    for y in x.sub_items
    if condition4(y)
]  # WRONG - too complex
```

**Fix:**
```python
result = []
for x in items:
    if condition1(x) and condition2(x) and condition3(x):
        for y in x.sub_items:
            if condition4(y):
                result.append(transform(x))
```

## 16. Not Using Path Objects

**Anti-pattern:**
```python
import os
path = os.path.join(dir_name, "file.txt")
if os.path.exists(path):
    with open(path) as f:
        ...
```

**Fix:**
```python
from pathlib import Path
path = Path(dir_name) / "file.txt"
if path.exists():
    with path.open() as f:
        ...
```

## 17. String Formatting with + or %

**Anti-pattern:**
```python
message = "Hello, " + name + "! You have " + str(count) + " messages."
message = "Hello, %s! You have %d messages." % (name, count)
```

**Fix:**
```python
message = f"Hello, {name}! You have {count} messages."
```

## 18. Not Using dataclasses

**Anti-pattern:**
```python
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def __repr__(self):
        return f"Point(x={self.x}, y={self.y})"
    
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y
```

**Fix:**
```python
from dataclasses import dataclass

@dataclass
class Point:
    x: float
    y: float
```

## 19. Lambda Abuse

**Anti-pattern:**
```python
process = lambda x: x.strip().lower().replace(" ", "_")[:20]  # WRONG
```

**Fix:**
```python
def process(x: str) -> str:
    """Clean and truncate string."""
    return x.strip().lower().replace(" ", "_")[:20]
```

## 20. Not Using Sets for Membership Testing

**Anti-pattern:**
```python
valid_codes = ["A1", "A2", "A3", ...]  # Long list
if code in valid_codes:  # O(n) lookup
    ...
```

**Fix:**
```python
valid_codes = {"A1", "A2", "A3", ...}  # Set
if code in valid_codes:  # O(1) lookup
    ...
```

## 21. Using Global Variables for State

**Anti-pattern:**
```python
current_user = None  # Global state

def login(username):
    global current_user  # WRONG - hard to test, track, and debug
    current_user = username

def get_current_user():
    return current_user
```

**Fix:**
```python
@dataclass
class AppState:
    current_user: str | None = None

# Pass state explicitly or use dependency injection
def login(state: AppState, username: str) -> None:
    state.current_user = username

# Or use a class
class UserSession:
    def __init__(self):
        self._current_user: str | None = None

    def login(self, username: str) -> None:
        self._current_user = username

    @property
    def current_user(self) -> str | None:
        return self._current_user
```

## 22. Star Imports

**Anti-pattern:**
```python
from os.path import *  # WRONG - pollutes namespace, unclear what's imported
from mymodule import *  # WRONG - can shadow builtins, hard to track

join("a", "b")  # Where does this come from?
```

**Fix:**
```python
from os.path import join, exists, dirname
# Or
import os.path

os.path.join("a", "b")  # Clear origin
```

**Why it's wrong:**
- Pollutes the namespace with unknown names
- Can shadow builtins or other imports
- Makes code hard to read and debug
- Breaks static analysis tools

## 23. Circular Imports

**Anti-pattern:**
```python
# module_a.py
from module_b import func_b  # WRONG - circular dependency

def func_a():
    return func_b() + 1

# module_b.py
from module_a import func_a  # ImportError!

def func_b():
    return func_a() + 1
```

**Fixes:**

1. **Import inside function** (quick fix):
```python
# module_a.py
def func_a():
    from module_b import func_b  # Deferred import
    return func_b() + 1
```

2. **Restructure modules** (proper fix):
```python
# shared.py - extract shared functionality
def base_func():
    return 1

# module_a.py
from shared import base_func

def func_a():
    return base_func() + 1

# module_b.py
from shared import base_func

def func_b():
    return base_func() + 2
```

3. **Use TYPE_CHECKING for type hints**:
```python
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from module_b import ClassB  # Only imported for type checking

def func_a(obj: ClassB) -> int:  # Works at runtime
    ...
```

## 24. Catching and Ignoring Exceptions

**Anti-pattern:**
```python
try:
    result = risky_operation()
except Exception:
    pass  # WRONG - silently hides all errors

try:
    data = fetch_data()
except:  # WRONG - catches KeyboardInterrupt, SystemExit
    return None
```

**Fix:**
```python
try:
    result = risky_operation()
except SpecificError as e:
    logger.warning("Operation failed, using default: {}", e)
    result = default_value

# If you must catch broadly, at least log
try:
    data = fetch_data()
except Exception:
    logger.exception("Unexpected error fetching data")
    raise  # Re-raise after logging
```

## 25. Hardcoding Paths

**Anti-pattern:**
```python
config = open("/home/user/project/config.json")  # WRONG - breaks on other machines
output = open("C:\\Users\\name\\output.txt")  # WRONG - Windows-specific
```

**Fix:**
```python
from pathlib import Path

# Relative to script/module
CONFIG_PATH = Path(__file__).parent / "config.json"

# Or use environment/config
from os import environ
DATA_DIR = Path(environ.get("DATA_DIR", "./data"))
output_path = DATA_DIR / "output.txt"
```

## Summary

Key principles to avoid anti-patterns:

1. Use built-in functions and standard library when possible
2. Leverage context managers for resource management
3. Use appropriate data structures (sets for membership, Counter for counting)
4. Keep code readable and idiomatic
5. Use modern Python features (f-strings, dataclasses, Path)
6. Avoid premature optimization
7. Write explicit, clear code over clever code
