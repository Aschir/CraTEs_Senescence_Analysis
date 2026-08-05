# GitHub Upload Guide

This folder is a clean analysis repository export. The original working folder
is still connected to `https://github.com/pulvercyril/crates`; do not push this
analysis directly to that remote.

## 1. Check The Export

From the terminal:

```bash
cd /Users/vinzent/craTEs/CraTEs_Senescence_Analysis
find . -type f -size +90M
git status -sb
```

The `find` command should print nothing for a normal GitHub upload. If it prints
a file, do not commit that file unless you use Git LFS.

## 2. Create A New GitHub Repository

Recommended repository name:

```text
CraTEs_Senescence_Analysis
```

Create an empty GitHub repository in the browser. Do not add a README, license,
or `.gitignore` in the browser because this export already has them.

## 3. Push The Existing Local Commit

This export has already been initialized as a local Git repository with an
initial commit.

After creating the empty GitHub repository, run:

```bash
cd /Users/vinzent/craTEs/CraTEs_Senescence_Analysis
git remote add origin https://github.com/<YOUR_GITHUB_USER>/CraTEs_Senescence_Analysis.git
git push -u origin main
```

Replace `<YOUR_GITHUB_USER>` with your GitHub username.

If you want the commit author to use a different email before pushing, run:

```bash
git config user.name "Vinzent Aschir"
git config user.email "<YOUR_GITHUB_EMAIL>"
git commit --amend --reset-author --no-edit
```

Then push as shown above.

## 4. Optional: Use GitHub CLI Later

`gh` is not installed on this machine right now. On macOS you can install and
authenticate it with:

```bash
brew install gh
gh auth login
```

Then a future repository can be created and pushed with:

```bash
gh repo create CraTEs_Senescence_Analysis --private --source . --remote origin --push
```

Use `--public` instead of `--private` only when the code, outputs, and data
access notes are ready to be public.

## 5. Attribution Rules

- Keep upstream craTEs source files in `R/`.
- Keep the craTEs MIT license in `LICENSES/craTEs-MIT.md`.
- Keep local changes and analysis-specific functions in
  `local_craTEs_extension/` or `scripts/`.
- Keep `scripts/install_local_crates.R` as the local packaging bridge that lets
  a fresh checkout install the vendored upstream code as `library(craTEs)`.
- If you modify upstream craTEs source directly, add a clear file header saying
  what was changed and by whom.

## 6. Result Files

Normal GitHub rejects files larger than 100 MB. The large
`target_gene_matrix_*.tsv` files are therefore ignored. If you need to publish
them, use Git LFS or archive them externally and link them from
`data/data_acces.md` or the README.
