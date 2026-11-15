# Publishing this repository to GitHub

A full walkthrough, written for a first repository. Delete this file once the
repository is live.

---

## 0. Before anything

**This repository now contains your melittin HDX-MS measurements**
(`data/hdx/melittin_observed.csv`). That was a deliberate choice. Two things
follow from it:

- Confirm with your supervisor, and with anyone who ran or co-owns those
  measurements, that publishing them ahead of the paper is agreed. The code is
  yours to publish; the data may have more than one stakeholder.
- If you change your mind, uncomment the two lines under "Experimental data" in
  `.gitignore` **before** your first commit. Once something is committed and
  pushed it stays in the history even if you delete the file later, and removing
  it properly means rewriting history.

Also check the licence. `LICENSE` is MIT, which is permissive and standard for
research software. Check KCL and funder policy before committing.

---

## 1. One-time setup

Check whether you already have git:

```bash
git --version
```

If not, install it. On macOS, `xcode-select --install` is the simplest route.

Tell git who you are. This is stamped on every commit and is public:

```bash
git config --global user.name "Alexandros Gkounta"
git config --global user.email "you@example.com"
```

Use an address you are happy to have public, or GitHub's privacy address
(Settings, Emails, "Keep my email addresses private") which looks like
`12345678+username@users.noreply.github.com`.

Set the default branch name to `main`:

```bash
git config --global init.defaultBranch main
```

---

## 2. Authentication

GitHub stopped accepting passwords for git operations in 2021. Pick one route.

### Option A: GitHub CLI (easiest)

```bash
brew install gh        # or see cli.github.com
gh auth login          # choose HTTPS, authenticate in the browser
```

This handles credentials for you and can create the remote repository from the
terminal. If you are unsure, use this one.

### Option B: SSH key

```bash
ssh-keygen -t ed25519 -C "you@example.com"      # press Enter at each prompt
cat ~/.ssh/id_ed25519.pub                        # copy the whole line
```

Paste it into GitHub, Settings, SSH and GPG keys, New SSH key. Test it:

```bash
ssh -T git@github.com
```

"Hi username! You've successfully authenticated" means it worked.

### Option C: Personal access token

GitHub, Settings, Developer settings, Personal access tokens, Fine-grained
tokens. Give it Contents read and write on the repository. Use the token in
place of a password when git prompts. Store it in your keychain so you are not
pasting it repeatedly.

---

## 3. Create the repository on GitHub

**Create it empty.** No README, no .gitignore, no licence, because this folder
already has all three and GitHub's versions would collide.

With the CLI, skip ahead; `gh repo create` in step 5 does it for you.

In the browser: github.com, New repository, name it `amp-md-pipeline`,
description "Reproducible MD processing and analysis for membrane-active
antimicrobial peptides, with HDX-MS forward modelling", Public, and create
without any initialisation files.

---

## 4. Commit locally

```bash
cd ~/Desktop/AMP_MD_processed/amp-md-pipeline
git init
git add .
```

**Now stop and read what you are about to commit:**

```bash
git status
```

Expect roughly 52 files. Check specifically that:

- No `.xtc`, `.tpr`, `.trr` or `.edr` files appear. Trajectories are large and
  do not belong in git.
- `data/hdx/melittin_observed.csv` appears **only if** you intend to publish it.
- No `__pycache__`, `.bak`, `.DS_Store` or `._` files appear.

If something is listed that should not be, add it to `.gitignore`, then:

```bash
git rm -r --cached .
git add .
git status
```

When the list is right:

```bash
git commit -m "Initial commit: AMP MD processing and analysis pipeline"
```

---

## 5. Push

With the CLI, this creates the remote and pushes in one step:

```bash
gh repo create amp-md-pipeline --public --source=. --remote=origin --push
```

Otherwise, having created the empty repository in step 3:

```bash
git remote add origin git@github.com:<username>/amp-md-pipeline.git   # SSH
# or
git remote add origin https://github.com/<username>/amp-md-pipeline.git   # HTTPS

git branch -M main
git push -u origin main
```

Reload the GitHub page. The README renders on the front page.

---

## 6. Day to day, from now on

```bash
git status                       # what has changed
git diff                         # the actual changes, line by line
git add py/amp_analysis.py       # stage one file
git add .                        # or everything
git commit -m "Fix contact cutoff in membrane analysis"
git push
```

Commit messages: say what changed and why, in the imperative. "Fix transposed
OBS indexing in HDX loader" is useful in a year. "update" is not.

Commit often and in small pieces. A commit per logical change is far easier to
read back than one enormous commit a month.

---

## 7. When something goes wrong

| Message | What it means | Fix |
| --- | --- | --- |
| `remote origin already exists` | You ran `git remote add` twice | `git remote set-url origin <url>` |
| `Updates were rejected because the remote contains work that you do not have` | You initialised the GitHub repository with a README | `git pull --rebase origin main` then push, or delete and recreate it empty |
| `Permission denied (publickey)` | SSH key not registered or not loaded | Re-run step 2B, check `ssh -T git@github.com` |
| `Support for password authentication was removed` | Using a password over HTTPS | Use a token or the CLI, step 2 |
| `fatal: not a git repository` | Wrong directory | `cd` to the repository folder |
| A large file is rejected | A trajectory slipped past `.gitignore` | Add the pattern, `git rm --cached <file>`, commit again |

**Committed something you should not have?** If you have not pushed yet,
`git reset --soft HEAD~1` undoes the commit and keeps your files. If you have
pushed, the file is in the history permanently until the history is rewritten,
which is a separate and more involved job. This is why step 4 says stop and read
`git status`.

---

## 8. After it is live

- [ ] Add topics on GitHub: `molecular-dynamics`, `gromacs`, `mdanalysis`,
      `hdx-ms`, `antimicrobial-peptides`, `computational-biophysics`.
- [ ] Replace `github.com/[USERNAME]` in your CV header and the Projects entry
      with the real URL.
- [ ] Add it to the Featured section of your LinkedIn profile.
- [ ] Consider a Zenodo DOI once the paper is submitted, so the code can be
      cited independently of the paper.

### Worth doing next, but not blocking

- A worked example on a public structure, so a visitor can run something without
  your trajectories.
- A few tests around `kint.py` and the descriptor functions. They are the most
  testable parts and the ones a reviewer would probe first.
