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
- **`ElementSymbol`** – AMR-specific search mode 

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

`FGRegex` allows users to define **regular-expression–like patterns** over compact Pfam architecture strings.  
It is primarily used to identify **order-sensitive domain architectures**, such as operon-like gene arrangements or specific domain compositions.

Example compact FG annotation of a protein:

`PF00001:::PF00002:::PF00003`


An FGRegex query can express constraints on **order**, **repetition**, and **optionality** of Pfam domains.

---

### 4.2 Query Format

Each FGRegex query is specified in a two-column, tab-separated format:



`protein_id<tab>FGRegex`
`query_1<tab>*, PF05658, *, PF03895`
`query_2<tab>(PF01234 | PF05678){1,5}, PF09012`


---

### 4.3 FGRegex Syntax Reference

| Syntax        | Meaning                                      | Example                      |
|---------------|----------------------------------------------|------------------------------|
| `*`           | Matches any number of arbitrary domains      | `*, PF00001`                 |
| `(A \| B)`    | Matches domain A **or** domain B              | `(PF00001 \| PF00002)`       |
| `{n,m}`       | Repeats the preceding domain *n* to *m* times | `PF00001{1,3}`               |
| `,`           | Domain separator (removed during parsing)    | `PF00001, PF00002`           |

---

### 4.4 Matching Rules and Mandatory-Domain Resolution

FGRegex matching is performed in two logically distinct stages:

1. **Pattern compilation and regex matching**
2. **Specificity-based conflict resolution using mandatory-domain counts**

The implementation of these stages is distributed across the protein-search layer and the graph-labeling layer.

---

#### Pattern Compilation and Matching

**Code location:**  
- `OperonData.convert_regex_fg_to_pattern()`  
- `single_cohort_API.compile_fg_query_pattern()`  
- `single_cohort_API.matches_fg_query()`

**Behavior:**

- Each FGRegex query is converted into a Python regular-expression pattern.
- During graph construction, compiled patterns are cached together with their metadata (including mandatory-domain count).
- A protein or graph node is considered a match if its **compact FG string**, with domain order preserved, satisfies the compiled regex.

> Note:  
> In the current implementation, regex compilation occurs **once per query during setup**, not globally at application startup.

---

#### Mandatory Domain Definition

**Code location:**  
- `OperonData.count_mandatory_domains()`  
- (protein-search anchoring only: `single_cohort_API.get_mandatory_domains()`)

Each FGRegex query is analyzed to determine its **mandatory-domain count**, which is used as a proxy for query specificity.

A Pfam domain contributes **+1** to the mandatory-domain count if and only if:

- It appears explicitly as a Pfam identifier (`PF\d{4,6}`),
- It is not a pure wildcard segment (`*`, `.`, or `.*`),
- It is not part of an alternation group (`|`),
- It is not optional due to a quantifier:
  - `?`
  - `*`
  - `{0,n}` (i.e., minimum repetition count is zero).

All other domains are treated as optional and **do not** increase specificity.

> Implementation note:  
> `OperonData.count_mandatory_domains()` returns an **integer count**, not a set, and does **not deduplicate** repeated mandatory segments.

---

#### Conflict Resolution Between Overlapping Matches

**Code location:**  
- `OperonData.get_matching_advanced_regex_queries()`  
- `OperonData.merge_with_prefix()`

Multiple FGRegex queries may match the same protein or graph node.  
To ensure deterministic and reproducible labeling, conflicts are resolved as follows:

1. For each matching FGRegex query, retrieve its precomputed mandatory-domain count.
2. Retain only the query (or queries) with the **maximum mandatory-domain count**.
3. If a single query remains, its label is assigned directly.
4. If multiple queries remain with equal maximal specificity:
   - Their labels are merged using a **longest-common-prefix–based strategy**.
   - The shared prefix is preserved, and diverging suffixes are combined using `/`.

This logic is executed during node labeling in the graph construction phase and guarantees consistent resolution of overlapping FGRegex matches.

---

#### Relation to Database-Level Search Optimization (Clarification)

Mandatory domains are also used earlier in the pipeline for **database-level candidate anchoring**:

- **Code location:** `single_cohort_API.search_protein_by_FG()`
- Mandatory Pfam domains are extracted via `get_mandatory_domains()`
- The rarest mandatory domain is selected to restrict the candidate protein set before regex evaluation

This optimization is **independent** of the conflict-resolution logic described above and only affects protein retrieval performance, not label assignment.

---

## 5. FGPerm: Permutation-Based (Order-Independent) Search

### 5.1 Motivation

`FGPerm` is designed for cases where **the presence of a set of Pfam domains matters**, but **their relative order does not**.  
This is common for modular proteins, mobile genetic elements, and recombined functional modules.

In contrast to `FGRegex`, `FGPerm` treats each query as an **unordered set of mandatory domains** and performs **set-containment matching**.

---

### 5.2 Query Format

Each FGPerm query is specified as:



`protein_id<tab>FGPerm`
`query_1<tab>PF01234,PF05678,PF09012`
`query_2<tab>PF00001,PF00002`

**Semantics:**

- Domains are separated by commas.
- Domain order is ignored.
- **All listed domains are mandatory**.
- No wildcard, alternation, or quantifier syntax is supported.

---

### 5.3 Matching Semantics

Let:

- Query domain set:  
  `Q = {PF00001, PF00003}`
- Protein FG annotation:  
  `PF99999:::PF00003:::PF00001:::PF88888`

The protein matches the query if and only if:


`Q ⊆ set(FG_annotation)`



No assumptions are made about:

- domain order,
- adjacency,
- or spacing.

**Canonical implementation:**

- `utils.has_all_mandatory_domains(mandatory_domains, protein_fg_string)`

This function provides the authoritative definition of FGPerm matching by converting the FG string into a set and testing set containment.

---

### 5.4 Database-Level Execution Strategy

To ensure scalability on large protein databases, FGPerm avoids full-database scans by using a **two-step candidate filtering strategy**.

This logic is implemented in:

- `single_cohort_API.search_proteins_by_mandatory_domains(...)`

---

#### Step 1: Candidate Anchoring (Selectivity Optimization)

**Code locations:**

- `single_cohort_API.search_proteins_by_mandatory_domains(...)`
- `single_cohort_API.count_prots_per_pfam(...)`

**Procedure:**

1. For each mandatory Pfam domain in the query:
   - Count how many proteins in the database contain that domain.
2. Identify the **rarest mandatory domain** (lowest protein count).
3. Retrieve **only proteins containing this domain** as initial candidates.

This ensures that the most selective constraint is applied first and dramatically reduces the search space.

---

#### Step 2: Set-Containment Filtering

**Code location:**

- `single_cohort_API.search_proteins_by_mandatory_domains(...)`
- `utils.has_all_mandatory_domains(...)`

**Procedure:**

- For each candidate protein:
  - Fetch its compact FG annotation string from the database.
  - Convert the FG string into a set of Pfam domains.
- Retain the protein if and only if: `mandatory_domain_set ⊆ protein_domain_set`


Only proteins satisfying **all** mandatory domains survive this step.

As a result, runtime scales with query selectivity rather than total database size.

---

### 5.5 Node-Level Matching (Graph Annotation)

For graph construction and node annotation, FGPerm operates directly on node-level FG annotations.

**Code location:**

- `OperonData.add_protein_ids()`

**Procedure:**

- The node’s FG annotation string is split into a set of Pfam domains.
- Each FGPerm query matches the node if: `mandatory_domain_set ⊆ node_domain_set`

- All matching queries are collected for subsequent conflict resolution.

---

### 5.6 Conflict Resolution and Specificity Ranking

Multiple FGPerm queries may match the same protein or graph node.  
To ensure deterministic and biologically meaningful labeling, conflicts are resolved using **mandatory-domain specificity**.

**Code locations:**

- `OperonData.get_matching_permutation_queries()`
- `OperonData.merge_with_prefix()`

**Resolution rules:**

1. For each matching query, compute its **mandatory-domain count**  
 (defined as the number of Pfam domains listed in the FGPerm query).
2. Retain only the query (or queries) with the **maximum mandatory-domain count**.
3. If a single query remains:
 - Assign its associated `protein_id`.
4. If multiple queries remain with equal maximal specificity:
 - Merge their labels deterministically using a **longest-common-prefix–based strategy**.

This guarantees that:

- More specific functional hypotheses take precedence over more permissive ones.
- Ambiguity is preserved only when specificity is genuinely indistinguishable.

---

### 5.7 Relation to FGRegex (Clarification)

- FGPerm uses **unordered set containment**.
- FGRegex uses **order-sensitive regex matching**.
- Both methods share the same **specificity-based conflict resolution strategy**, but differ in:
  - how mandatory domains are defined,
  - how matches are computed,
  - and how database-level candidate reduction is performed.

---

## 6. AMR Search: Antimicrobial Resistance Annotation (ElementSymbol)

In addition to domain-based functional group searches, the pipeline supports an **AMR-specific search and annotation mode**, exposed as the search type **`ElementSymbol`**.  
This mode integrates antimicrobial resistance (AMR) information into MAGraph as an **orthogonal annotation layer**, rather than as a competing domain query language.

AMR search is designed to answer a different class of biological questions:

> *Which genomic neighborhoods contain known antimicrobial resistance genes, and how are they embedded in functional and genomic context?*

---

### 6.1 Conceptual Role of AMR Search

Unlike `FG`, `FGRegex`, or `FGPerm`, AMR search:

- is **not pattern-based**,
- does **not operate on Pfam domain architectures**,
- does **not participate in query specificity ranking or conflict resolution**.

Instead, AMR search provides **direct gene-level annotations** derived from external AMR detection pipelines (e.g. CARD, AMRFinderPlus), which are then projected onto MAGraph nodes.

Conceptually:

| Aspect | FG / FGRegex / FGPerm | AMR (ElementSymbol) |
|------|------------------------|---------------------|
| Input | User-defined queries | External AMR hits |
| Semantics | Domain-based matching | Resistance gene annotation |
| Matching logic | Regex / set containment | Exact hit mapping |
| Conflict resolution | Yes (mandatory-domain specificity) | No |
| Role in pipeline | Drives graph selection | Enriches node annotation |

---

### 6.2 Search Type Declaration

AMR search is activated by declaring the search type as **`ElementSymbol`** in the input file header.

Example (`functional_groups.txt`):

`protein_id ElementSymbol`
`blaTEM`
`tetM`
`ermB`


- Each entry represents an **AMR element or resistance symbol**.
- The interpretation of symbols depends on the upstream AMR annotation source.
- No domain syntax, wildcards, or regex constructs are supported.

---

### 6.3 Data Source and Matching Semantics

**Matching semantics:**

- A protein (or gene cluster) matches an AMR query if it has been annotated with the corresponding AMR symbol.
- Matching is **binary**: a gene either carries an AMR annotation or it does not.
- Multiple AMR annotations may coexist on the same gene cluster.

**Key characteristics:**

- No mandatory-domain concept
- No query competition
- No specificity ranking

AMR annotations are therefore **non-exclusive** and purely descriptive.

---

### 6.4 Implementation Details

**Core code locations:**

- `operon_data.OperonData.add_amr_ids()`
- `operon_data.OperonData.amr_hit_dic`

**Integration point:**

- AMR annotations are applied during the **node annotation phase**, after MAGraph construction.

**Typical flow:**

1. External AMR detection produces mappings of the form: `protein_id → AMR_symbol`
2. These mappings are loaded into `OperonData`.
3. During `add_protein_ids()`:
- AMR symbols are attached to the corresponding graph nodes (`gc_id`).
4. Nodes may accumulate:
- FG / FGRegex / FGPerm labels **and**
- one or more AMR labels.

---

### 6.5 Interaction with MAGraph Profiling

AMR annotations do **not** influence:

- protein search results,
- gene cluster selection,
- conflict resolution between FG queries.

However, they are fully available for:

- node coloring and visualization,
- subgraph interpretation,
- downstream filtering (e.g. “subgraphs containing at least one AMR gene”),
- resistance-focused comparative analyses.

AMR therefore acts as a **semantic overlay** on top of the domain-driven MAGraph structure.

---

### 6.6 Relation to Other Search Modes

AMR search complements, rather than replaces, functional group searches:

- FG / FGRegex / FGPerm determine *where* graphs are built.
- AMR annotations help explain *what functional risks or traits* are present within those graphs.

This separation preserves clarity of semantics while enabling rich, multi-layered genomic neighborhood analysis.

---


## 7. Integration into the MAGraph Pipeline

This section describes how **domain-based searches** (`FG`, `FGRegex`, `FGPerm`) and **AMR annotation** (`ElementSymbol`) are integrated into the MAGraph construction, annotation, and profiling workflow, with clarifications aligned to the actual implementation.

---

### 7.1 Cohort Mapping (Query-to-Gene-Cluster Resolution)

All **domain-based search modes** (`FG`, `FGRegex`, `FGPerm`) ultimately produce a unified mapping of the form:

`{ query → [gc_id, ...] }`


**Code locations:**

- `parser.load_domain_search_result_and_search_type()`
- `utils.cohort_domain_search`
- `single_cohort_API.get_FG2prot_id(...)`
- `single_cohort_API.get_FG2prot_id_regex(...)`
- `single_cohort_API.search_proteins_by_mandatory_domains(...)`

**Details:**

- Each query is associated with a set of matching **protein IDs**, which are then mapped to **gene cluster IDs (`gc_id`)**.
- The mapping is normalized across all domain-based search modes so downstream steps remain agnostic to the original query type.
- Internally, both forward mappings (`query → gc_ids`) and reverse mappings (`gc_id → query`) are constructed.

This unified mapping serves as the foundation for:
- node labeling and coloring,
- query coverage and hit statistics,
- subsequent graph construction and filtering.

> **AMR clarification:**  
> `ElementSymbol` (AMR) does **not** produce `{query → gc_id}` mappings.  
> Instead, AMR annotations are applied directly at the **node annotation stage** (see §7.3).

---

### 7.2 Graph Construction (MAGraph Assembly)

Matched `gc_id`s are expanded into **gene neighborhood graphs (MAGraphs)**.

**Code locations:**

- `retrieve_graph_extraction_information.py`
- `graph_classes.py`
- `runner.py`

**Graph model:**

- **Nodes** represent gene clusters (`gc_id`).
- **Edges** represent genomic adjacency relationships.
- **Edge attributes** encode:
  - co-occurrence support across samples,
  - strand consistency,
  - optional distance or orientation metadata.

Graph construction is **query-driven**: only gene clusters connected to at least one matched `gc_id` (from FG / FGRegex / FGPerm) are expanded, ensuring scalability.

> AMR annotations do not influence which graphs are constructed.

---

### 7.3 Node Annotation, Conflict Resolution, and AMR Integration

After graph construction, nodes are annotated with functional and resistance-related labels.

**Code locations:**

- `operon_data.OperonData.add_protein_ids()`
- `operon_data.OperonData.add_amr_ids()`

**Behavior:**

- Each node may accumulate:
  - **domain-based labels** (`FG`, `FGRegex`, `FGPerm`), and
  - **AMR annotations** (`ElementSymbol`).
- For `FGRegex` and `FGPerm`, **node-level conflict resolution** is applied:
  - queries are ranked by **mandatory-domain specificity**,
  - only maximally specific queries are retained,
  - ties are resolved via longest-common-prefix–based label merging.
- **AMR annotations are non-competitive**:
  - multiple AMR labels may coexist on the same node,
  - no specificity ranking or conflict resolution is applied.

This design ensures:
- deterministic functional labeling,
- orthogonal integration of resistance information.

---

### 7.4 Subgraph Profiling and Ranking

Constructed MAGraphs (or their connected components / subgraphs) are profiled and ranked.

**Code locations:**

- `graph_classes.py`
- `operon_data.py`

**Ranking criteria include:**

- number of matched **domain queries** per subgraph,
- cumulative or median abundance (e.g., RPKM) across matched queries,
- prevalence across samples or cohorts.

**Optional annotations:**

- MSP clustering or GTDB-Tk–derived taxonomy,
- functional group summaries aggregated from node labels,
- AMR label presence summaries (descriptive, not ranking-driving).

These metrics enable prioritization of biologically meaningful gene neighborhoods while retaining resistance context.

---

### 7.5 Execution Order and Separation of Concerns (Clarification)

The MAGraph pipeline enforces a clear separation between stages:

1. **Protein-level search**  
   (`FG`, `FGRegex`, `FGPerm`)
2. **Gene-cluster mapping**
3. **Graph construction**
4. **Node annotation**  
   - domain-based labeling + conflict resolution  
   - AMR annotation (ElementSymbol)
5. **Subgraph profiling and optional constraint filtering**

This separation ensures that:
- search semantics are isolated from graph topology,
- AMR annotations do not bias domain-driven search results,
- new query or annotation modes can be added without modifying graph logic,
- profiling strategies can evolve independently of matching algorithms.

---

## 8. Summary

This pipeline provides a flexible, extensible framework for functional and resistance-aware exploration of MAGraph data.

In particular:

- **FGRegex** supports expressive, order-aware Pfam domain patterns.
- **FGPerm** enables robust, order-independent domain set queries.
- **AMR (ElementSymbol)** provides orthogonal antimicrobial resistance annotation at the node level.

Together, these components allow users to move seamlessly from:

> **protein domain logic**  
> → **cohort-scale genomic neighborhood structure**  
> → **functional and AMR-aware biological interpretation**.



