# School Study Materials Generator

Generate Kurzusammenfassungen (summaries), Lernziele (learning objectives), and Kontrollfragen (control questions) from school scripts.

## Workflow Overview

```
PDF (Original Script)
     │
     ▼
 Markdown (_EXT.md)
     │
     ├──► Lernziele.md (Learning Objectives)
     │        │
     │        ▼
     ├──► Spickzettel/Kurzusammenfassung.md (Cheat Sheet)
     │
     └──► Kontrollfragen.md (Practice Questions)
```

## Step 1: PDF to Markdown Conversion

### Source Materials
- Original school scripts in `orig/` folder
- Extract images to `images/` folder organized by document

### Naming Convention
| Suffix | Purpose |
|--------|---------|
| `_ORIG.pdf` | Original school document |
| `_EXT.md` | Extended version with explanations |
| `_loesung.md` | Version with solutions/answers |

### Markdown Structure

```markdown
# Hauptthema

## Lernziele

**Lernziel 1:** [Beschreibung]
**Lernziel 2:** [Beschreibung]

---

## Thema 1

### Unterthema 1.1

Inhalt...

---

## Zusammenfassung

[Tabellen und Übersichten]
```

## Step 2: Pandoc Formatting Guidelines

### Makefile Template

```makefile
# Makefile for converting Markdown to PDF
# Requires: pandoc, xelatex (MacTeX/TinyTeX)

MD_FILES := $(wildcard *.md)
PDF_FILES := $(MD_FILES:.md=.pdf)

PANDOC_OPTS := --pdf-engine=xelatex \
               -V geometry:margin=2.5cm \
               -V fontsize=11pt \
               -V lang=de \
               -V mainfont="Helvetica" \
               --toc \
               --toc-depth=3 \
               -V header-includes='\usepackage{array}\renewcommand{\arraystretch}{1.3}'

all: $(PDF_FILES)

%.pdf: %.md
	pandoc $(PANDOC_OPTS) -o $@ $<

clean:
	rm -f $(PDF_FILES)

.PHONY: all clean
```

### Chemical Formulas (Subscripts)

Use pandoc subscript notation with tildes:
- `H~2~O` renders as H₂O
- `CO~2~` renders as CO₂
- `C~6~H~12~O~6~` renders as C₆H₁₂O₆

### Mathematical Formulas

Use LaTeX math notation:
- Inline: `$\alpha = \beta$`
- Display: `$$\alpha = \beta$$`

### Table Formatting

**DO:**
- Keep headers short (1-3 words)
- Use k-notation for large numbers (300k, 1.9 Mio)
- Avoid bold text in table cells
- Use consistent separator dashes

**DON'T:**
- No bold `**text**` inside tables
- No long header descriptions like "Wann (vor ... Jahren)"
- No Swiss number formatting with apostrophes (300'000)

**Good Example:**
```markdown
| Art | Wann | Wo | Besonderheiten |
|-----|------|-----|----------------|
| Homo erectus | 1.9 Mio - 110k | Afrika, Asien | Erste Art ausserhalb Afrikas |
```

**Bad Example:**
```markdown
| **Art** | **Wann (vor ... Jahren)** | **Wo gefunden** |
|---------|---------------------------|-----------------|
| **Homo erectus** | 1'900'000 - 110'000 | Afrika und Asien |
```

### Deep Nesting (Lists)

For documents with deeply nested lists, add header.tex:

```latex
% header.tex
\usepackage{enumitem}
\setlistdepth{9}
\setlist[itemize,1]{label=$\bullet$}
\setlist[itemize,2]{label=$\circ$}
\setlist[itemize,3]{label=$\diamond$}
\setlist[itemize,4]{label=$\ast$}
\renewlist{itemize}{itemize}{9}
```

And add to Makefile: `-H header.tex`

### Horizontal Rules

Use `---` to create visual separators between sections.

### ASCII Art Diagrams

Simple ASCII diagrams work well for concepts:

```markdown
```
    Windverbreitung              Tierverbreitung
         🌬️                          🐰
          |                           |
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```
```

## Step 3: Lernziele (Learning Objectives)

### Structure

Extract learning objectives from:
1. Teacher-provided Lernziele documents
2. Section headers in scripts
3. Implied objectives from content

### Format

```markdown
# Lernziele [Thema]

## Bereich 1

---

### Lernziel 1: [Kurztitel]

**Frage:** Ich kann [Kompetenz beschreiben].

**Antwort:**

[Ausführliche Erklärung mit Beispielen]

---

### Lernziel 2: [Kurztitel]

...
```

### Synchronization

IMPORTANT: Always verify that:
1. All Lernziele from the official list are covered
2. Content from all relevant scripts is included
3. Examples and details match the script content

## Step 4: Spickzettel/Kurzusammenfassung (Cheat Sheet)

### Purpose

Compact reference document with key facts for quick review before exams.

### Structure

```markdown
# Spickzettel: [Thema]

## 1. [Hauptbereich]

### Wichtige Begriffe

| Begriff | Definition |
|---------|------------|
| Term 1 | Kurze Definition |
| Term 2 | Kurze Definition |

### Wichtige Fakten

- Fakt 1
- Fakt 2
- Fakt 3

### Formeln/Tabellen

[Kompakte Übersichten]

---

## 2. [Nächster Bereich]

...

---

## Zusammenfassung

[Finale Übersichtstabelle mit allen wichtigen Punkten]
```

### Guidelines

- Maximum 2-3 pages per major topic
- Use tables for comparisons
- Include all required terminology
- Reference Lernziel numbers where relevant

## Step 5: Kontrollfragen (Practice Questions)

### Purpose

Q&A format for oral practice and self-testing.

### Structure

```markdown
# Kontrollfragen zu [Thema]

Frage-Antwort Format zum mündlichen Abfragen.

---

## Bereich 1 (LZ 1-5)

---

**F1:** [Frage]

**A:** [Antwort]

---

**F2:** [Frage]

**A:** [Antwort]

---

## Bereich 2 (LZ 6-10)

...

---

## Schnelle Wissensfragen

---

**F50:** [Kurze Frage]?

**A:** [Kurze Antwort]

---

## Zusammenfassung: Die wichtigsten Fakten

| Thema | Merksatz |
|-------|----------|
| Konzept 1 | Einfache Zusammenfassung |
| Konzept 2 | Einfache Zusammenfassung |
```

### Question Types

1. **Definition Questions**: Was ist [Begriff]?
2. **Comparison Questions**: Was ist der Unterschied zwischen X und Y?
3. **Explanation Questions**: Erkläre/Beschreibe [Konzept].
4. **Application Questions**: Wie würde sich X ändern wenn Y?
5. **List Questions**: Nenne drei Beispiele für [Kategorie].
6. **Quick Fire**: Short answer questions for rapid review

### Mapping to Lernziele

- Reference Lernziel numbers in section headers: `## Bereich (LZ 1-5)`
- Ensure every Lernziel has at least 2-3 related questions
- Include questions testing both recall and understanding

## Subject-Specific Notes

### Biologie

- Use proper biological terminology in German
- Include diagrams for processes (photosynthesis, evolution)
- Chemical formulas with subscript notation

### Physik

- Mathematical formulas with LaTeX
- ASCII diagrams for optics (light rays, shadows)
- Greek letters: $\alpha$, $\beta$, etc.

### Mathematik

- Step-by-step solutions (Rechenweg)
- Multiple difficulty levels
- Both DOCX and PDF output for flexibility

### Deutsch

- Grammar tables (Nomen, Verben)
- Example sentences
- Correction exercises (*_korrekturen.md)

## Quality Checklist

Before finalizing materials:

- [ ] All Lernziele from official list are covered
- [ ] Chemical/math formulas render correctly
- [ ] Tables are properly formatted (no bold, short headers)
- [ ] Horizontal rules separate sections
- [ ] PDF generates without errors
- [ ] Content matches original scripts
- [ ] No duplicate information
- [ ] Appropriate depth for grade level

## Directory Structure Example

```
Subject/
├── CLAUDE.md                  # Project documentation
├── Makefile                   # Build configuration
├── header.tex                 # Optional LaTeX customization
├── orig/                      # Original PDFs
│   ├── Script1_ORIG.pdf
│   └── Script2_ORIG.pdf
├── images/                    # Extracted images
│   ├── script1/
│   └── script2/
├── 01_Topic1_EXT.md          # Extended script
├── 01_Topic1_EXT.pdf
├── 02_Lernziele.md           # Learning objectives
├── 02_Lernziele.pdf
├── 03_Spickzettel.md         # Cheat sheet
├── 03_Spickzettel.pdf
├── 04_Kontrollfragen.md      # Practice questions
└── 04_Kontrollfragen.pdf
```

## Common Commands

```bash
# Build all PDFs
make all

# Build specific file
make 03_Spickzettel.pdf

# List available files
make list

# Clean generated files
make clean

# Preview PDF (macOS)
open 03_Spickzettel.pdf
```
