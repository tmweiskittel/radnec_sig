# radnec_sig

Reanalysis of the GBMatch glioblastoma methylation dataset to evaluate methylation patterns associated with radiation necrosis.

## Setup on Google Cloud VM

The analysis is run from the persistent disk mounted at `/home/jupyter`.

### Clone the repository

```bash
cd /home/jupyter

git clone https://github.com/tmweiskittel/radnec_sig.git
cd radnec_sig
```

Confirm the repository and working directory:

```bash
git status
pwd
ls -lah
```

### Download the Klughammer data

Create a local directory for the methylation data:

```bash
mkdir -p data/klughammer
```

Copy the data from Google Cloud Storage:

```bash
gcloud storage cp \
    --recursive \
    gs://weiskittel-projects1/radnecrosis/klughammer_data/* \
    data/klughammer/
```

Check the downloaded files and disk usage:

```bash
du -sh data/klughammer

find data/klughammer \
    -maxdepth 2 \
    -type f \
    -printf '%p\t%k KB\n' \
    | sort
```

The `data/` directory is excluded from version control through `.gitignore`.

## Analysis questions

Analyses will be restricted to samples labeled as **GB classic** where appropriate.

### 1. Radiation necrosis versus recurrent tumor

Compare samples diagnosed with radiation necrosis against recurrent samples without annotated radiation necrosis.

### 2. Type of necrosis in recurrent samples

Compare recurrent samples classified as:

- radiation necrosis
- pseudopalisading necrosis
- neither

### 3. Quantitative necrosis phenotype

Within samples annotated as having radiation necrosis, evaluate necrosis as a continuous variable.

As an initial step, assess correlations among the available necrosis-related histopathologic and imaging variables.
