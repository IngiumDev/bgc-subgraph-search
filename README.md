# Domain Search & MAGraph Profiling Pipeline

This repository implements a modular pipeline for **domain-based search**, **cohort-level gene cluster retrieval**, and **graph-based genomic neighborhood analysis** on MAGraph data. It supports multiple search paradigms ranging from exact functional-group lookup to advanced regex- and permutation-based domain queries, and integrates these seamlessly into downstream graph construction and profiling.

This README provides a **user-facing overview** of the pipeline, with particular focus on the **FGRegex (regular expression)** and **FGPerm (permutation) search modes**.

---

For specific topics, refer to:
* [Quick Start Guide](docs/quickstart.md)
* [Advanced Search (FGRegex/FGPerm) Examples](docs/search_examples.md)
* [Docs](docs/func_docs.md)

---

## 1. Pipeline Overview

At a high level, the pipeline consists of four stages:

1. **Domain Search**  
   Identify proteins / gene clusters that match user-defined functional group (FG) queries.

2. **Cohort-Level Mapping**  
   Map matched proteins to genomic context IDs (`gc_id`) across a cohort.

3. **MAGraph Construction**  
   Build gene neighborhood graphs based on physical adjacency and co-occurrence.

4. **Subgraph Profiling & Filtering**  
   Quantify abundance / prevalence and apply optional functional constraints.

Each stage is agnostic to the specific search strategy used, as long as it is declared consistently in the input.

---

## 2. Supported Domain Search Types

The domain search is controlled by a **search type** specified directly in the input file **`\code\constants\search_types.py`**. The following modes are supported:

- **`FG`** – Exact functional group match  
- **`FGRegex`** – Regex-based Pfam architecture match  
- **`FGPerm`** – Permutation (order-independent) functional group match  
- **`ElementSymbol`** – AMR-specific search mode (not covered here)

The search type determines how queries are interpreted and how matching is performed at both database and graph levels.

---

## 3. Input File Format (`functional_groups.txt`)

The domain search input file must be **tab-separated** and self-describing.

### 3.1 Header (mandatory)

The first line specifies the search type in the second column:

`protein_id FGPerm`

Valid values are: `FG`, `FGRegex`, `FGPerm`, `ElementSymbol`.

If the search type cannot be determined or is invalid, the pipeline will stop with an explicit error.

---

## 4. FGRegex: Regex-Based Domain Architecture Search

### 4.1 Concept

`FGRegex` allows users to define **regular-expression–like patterns** over compact Pfam architecture strings. This is useful for identifying operon-like structures or specific domain arrangements.

Example FG annotation of a protein:

`PF00001:::PF00002:::PF00003`

A regex query may express constraints on order, repetition, or optional domains.

### 4.2 Query Format

`protein_id FGRegex`
`queryA PF00001.*PF00003`
`queryB PF00567(PF00901)?`

### 4.3 Matching Rules

- Queries are compiled into regex patterns once.
- A protein matches if its FG string satisfies the regex.
- If multiple regex queries match the same protein or node:
  - The query with the **largest number of mandatory domains** is preferred.
  - Ties are resolved by merging query names using a common-prefix strategy.

---

## 5. FGPerm: Permutation-Based (Order-Independent) Search

### 5.1 Motivation

`FGPerm` is designed for cases where **the presence of a set of domains matters**, but **their order does not**. This is common for modular proteins, mobile genetic elements, and recombined functional modules.

Unlike `FGRegex`, `FGPerm` treats queries as **sets of mandatory domains**.

---

### 5.2 Query Format

`protein_id FGPerm`
`queryA PF00001,PF00003`
`queryB PF00567,PF00901,PF01234`

Semantics:
- Domains are separated by commas.
- Order is ignored.
- All listed domains are mandatory.

---

### 5.3 Matching Semantics

Let:

- Query set: `Q = {PF00001, PF00003}`
- Protein FG annotation:  `PF99999:::PF00003:::PF00001:::PF88888`
The protein matches the query if:
`Q ⊆ set(FG_annotation)`


No assumptions are made about order, adjacency, or spacing.

---

### 5.4 Database-Level Execution

To ensure scalability, `FGPerm` uses a **two-step strategy**:

1. **Candidate anchoring**  
   - Count how many proteins contain each mandatory domain.
   - Select the rarest domain.
   - Retrieve proteins containing that domain only.

2. **Set containment filtering**  
   - Fetch compact FG strings for candidate proteins.
   - Retain only those containing *all* mandatory domains.

This avoids full-database scans and keeps runtime manageable.

---

### 5.5 Conflict Resolution

If multiple `FGPerm` queries match the same protein or graph node:

- The query with the **largest number of mandatory domains** is selected.
- If still ambiguous, query names are merged deterministically.

---

## 6. Integration into the MAGraph Pipeline

### 6.1 Cohort Mapping

All search modes (`FG`, `FGRegex`, `FGPerm`) produce a unified mapping:

`{ query → [gc_id, ...] }`


This mapping is used for:
- node coloring
- query coverage statistics
- downstream graph filtering

---

### 6.2 Graph Construction

Matched `gc_id`s are expanded into **gene neighborhood graphs** (MAGraphs), where:
- nodes represent gene clusters
- edges represent genomic adjacency
- edge weights encode support and strand consistency

---

### 6.3 Subgraph Profiling

Subgraphs are filtered and ranked based on:
- number of matched queries
- abundance (median RPKM across queries)
- prevalence across samples

Optional extensions include:
- MSP / GTDB-Tk taxonomy annotation
- functional group summaries

---

### 6.4 Pfam-Constrained Subgraph Selection (Optional)

An optional profiling mode allows enforcing **required FGRegex / FGPerm constraints** at the subgraph level:

> A subgraph is retained only if **at least one node** matches a required domain query.

This enables detection of genomic neighborhoods that contain specific functional signatures.

---

## 7. Current Limitations / TODO

- [ ] Formal CLI documentation (arguments, examples)
- [ ] Unified visualization examples
- [ ] Prevalence-based graph scoring (currently partial)
- [ ] Performance benchmarks on large cohorts
- [ ] More worked examples for FGRegex vs FGPerm

(Feel free to ask or leave these sections blank for now.)

---

## 8. Summary

This pipeline provides a flexible, extensible framework for functional-group-driven exploration of MAGraph data. In particular:

- **FGRegex** supports expressive, order-aware domain patterns.
- **FGPerm** enables robust, order-independent domain set queries.
- Both integrate consistently into cohort search, graph construction, and subgraph profiling.

Together, these modes allow users to move seamlessly from **protein domain logic** to **cohort-scale genomic neighborhood analysis**.


