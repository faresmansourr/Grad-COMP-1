import pandas as pd
import requests
import time
import os
from tqdm import tqdm
import pickle

input_file = "kegg_output/taxon_to_ko.csv"
output_file = "kegg_output/taxon_kegg_pathways.csv"
pathway_name_file = "kegg_output/pathway_names.csv"
cache_file = "kegg_output/ko_to_pathways.pkl"

# === Step 1: Read KO mappings ===
df = pd.read_csv(input_file)
if "Taxon" not in df or "KO" not in df:
    raise ValueError("CSV must contain columns: Taxon, KO")

# === Step 2: Define KEGG functions ===
def get_pathways_for_ko(ko):
    url = f"http://rest.kegg.jp/link/pathway/ko:{ko}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            return [line.split("\t")[1].replace("path:", "") for line in r.text.strip().split("\n")]
    except:
        pass
    return []

def get_pathway_name(pathway_id):
    url = f"http://rest.kegg.jp/get/{pathway_id}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            for line in r.text.splitlines():
                if line.startswith("NAME"):
                    return line.replace("NAME", "").strip()
    except:
        pass
    return "Unknown"

# === Step 3: Map KO to Pathways (with caching) ===
ko_list = sorted(df["KO"].unique())
ko_to_pathways = {}

if os.path.exists(cache_file):
    with open(cache_file, "rb") as f:
        ko_to_pathways = pickle.load(f)
    print(f"✅ Loaded {len(ko_to_pathways)} cached KO-pathway mappings.")

print("🔗 Mapping KO terms to KEGG pathways...")
for i, ko in enumerate(tqdm(ko_list)):
    if ko in ko_to_pathways:
        continue
    try:
        ko_to_pathways[ko] = get_pathways_for_ko(ko)
    except Exception as e:
        print(f"[!] Error fetching {ko}: {e}")
        ko_to_pathways[ko] = []
    if i % 20 == 0:
        with open(cache_file, "wb") as f:
            pickle.dump(ko_to_pathways, f)
    time.sleep(0.5)

with open(cache_file, "wb") as f:
    pickle.dump(ko_to_pathways, f)
print("✅ All KO → Pathway mappings cached.")

# === Step 4: Fetch pathway names with progress + resume support ===
print("🧠 Fetching pathway names...")
pathway_names = {}
pathway_ids = sorted({pid for paths in ko_to_pathways.values() for pid in paths})

# Resume if file already exists
if os.path.exists(pathway_name_file):
    try:
        saved_df = pd.read_csv(pathway_name_file)
        pathway_names = dict(zip(saved_df["Pathway_ID"], saved_df["Pathway_Name"]))
        print(f"✅ Resuming from {len(pathway_names)} existing pathway names.")
    except:
        print("⚠️ Could not resume from saved file.")

for i, pid in enumerate(tqdm(pathway_ids)):
    if pid in pathway_names:
        continue
    try:
        pathway_names[pid] = get_pathway_name(pid)
    except Exception as e:
        print(f"[!] Error fetching {pid}: {e}")
        pathway_names[pid] = "Unknown"
    if i % 20 == 0:
        pd.DataFrame.from_dict(pathway_names, orient="index", columns=["Pathway_Name"])\
          .reset_index().rename(columns={"index": "Pathway_ID"})\
          .to_csv(pathway_name_file, index=False)
    time.sleep(0.3)

# Final save
pd.DataFrame.from_dict(pathway_names, orient="index", columns=["Pathway_Name"])\
    .reset_index().rename(columns={"index": "Pathway_ID"})\
    .to_csv(pathway_name_file, index=False)
print("✅ Pathway names saved.")

# === Step 5: Build final output ===
print("📦 Saving final Taxon-Pathway table...")
final = []
for _, row in df.iterrows():
    for pid in ko_to_pathways.get(row["KO"], []):
        final.append({
            "Taxon": row["Taxon"],
            "KO": row["KO"],
            "Pathway_ID": pid,
            "Pathway_Name": pathway_names.get(pid, "Unknown")
        })

pd.DataFrame(final).drop_duplicates().to_csv(output_file, index=False)
print("🎉 Done! Pathway output saved to:", output_file)
