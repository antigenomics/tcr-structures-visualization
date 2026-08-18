# tcr-structures-visualization

Generates the TCR–pMHC contact maps served by [VDJdb](https://vdjdb.com) on its structure pages.

For each modelled or native structure the pipeline produces a 2D contact map — the CDR3α, CDR3β and
peptide residues laid out by PCA projection of their Cα coordinates, with lines drawn between
residues in contact — plus the per-residue coordinates and the contact list behind it.

## Output

Per structure hash, written to the `vdjdb-db/structure/` tree and published to
[`isalgo/vdjdb_structure_models`](https://huggingface.co/datasets/isalgo/vdjdb_structure_models):

| File | What |
|---|---|
| `<hash>.svg` / `<hash>.html` | the contact map (matplotlib SVG, embedded in a bare HTML wrapper) |
| `<hash>_simplified.svg` / `.html` | reduced variant, used for stacked overlay layers |
| `<hash>_aa_coordinates.tsv` | `Chain, Residue, ResNum, X, Y, Z` for every residue |
| `<hash>_contacts_aa.txt` | `chain_from, aa_from, res_num_from, aa_to, res_num_to, chain_to` |

VDJdb fetches these with [`tools/sync_structure_files.py`](https://github.com/antigenomics/vdjdb-web/blob/master/tools/sync_structure_files.py)
and injects the HTML into the structure viewer.

## Layout

- `produce_plots_pipline/` — the pipeline. `run_pipeline.py` orchestrates; `plotting.py` draws the
  map and writes the sidecar files; `coordinates.py` and `alignment.py` are the supporting steps.
- `anastasia_notebooks/` — contact-analysis notebooks comparing contact-calling methods
  (5 Å cutoff vs GetContacts) across relaxed, Rosetta and TCRmodel structures.
- `utils_scripts_and_notebooks/` — assorted analysis and batch-run helpers.

## Authors

Pipeline and plotting code — **Daniil Luppov** ([@LuppovDaniil](https://github.com/LuppovDaniil)),
PhD student on the VDJdb project.

Contact-analysis notebooks (`anastasia_notebooks/`) — **Anastasiia Alexandrova**
([@nastyaleksa04](https://github.com/nastyaleksa04)).

Per-commit authorship is preserved in the git history, which was imported intact from the group's
internal GitLab.

## Known limitation

The emitted SVG carries only matplotlib's automatic element ids (`line2d_N`, `text_N`), so nothing
in the markup identifies which residue or which contact an element represents. That is what stops
the VDJdb viewer from highlighting a residue or a contact line on hover. The fix is to set a `gid`
on each artist as it is drawn — matplotlib writes an artist's `gid` straight out as the SVG `id`
attribute — at the two draw sites in `plotting.py` (backbone segments, and the Cα–Cα contact lines).
