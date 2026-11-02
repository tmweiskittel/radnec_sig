# Work from the persistent attached disk
cd /home/jupyter

# Clone the repository
git clone https://github.com/tmweiskittel/radnec_sig.git
cd radnec_sig

# Confirm repository
git status
pwd
ls -lah

# Create local data directory
mkdir -p data/klughammer

# Download Klughammer data from Google Cloud Storage
gcloud storage cp \
    --recursive \
    gs://weiskittel-projects1/radnecrosis/klughammer_data/* \
    data/klughammer/

# Check downloaded files and disk usage
du -sh data/klughammer
find data/klughammer -maxdepth 2 -type f -printf '%p\t%k KB\n' | sort


# radnec_sig
- limit to GB classic labeled samples 

1) Radionecrosis as diagnosis versus recurrence with no radionecrosis annotated on recurrence sample

2) Radionecrosis present vs pseudopalisading vs none in recurrence sample

3) Necrosis as continuous variable in samples labeled radionecrosis present
- check correlation of necrosis related variables. Anticipating they are all the same  
