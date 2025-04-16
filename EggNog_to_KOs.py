import os
import pandas as pd
import requests
import time

# === CONFIG ===
renamed_folder = "renamed_annotations"
output_folder = "kegg_output"

os.makedirs(output_folder, exist_ok=True)

# === 1. Read .tsv files with taxon names ===
tsv_files = sorted([f for f in os.listdir(renamed_folder) if f.endswith(".tsv")])
renamed_paths = [(os.path.splitext(f)[0], os.path.join(renamed_folder, f)) for f in tsv_files]

# === 2. Extract KO terms ===
records = []
skipped = []

print("🔍 Parsing renamed files and extracting KO terms...")
for taxon, file_path in renamed_paths:
    try:
        df = pd.read_csv(file_path, sep="\t", comment="#", low_memory=False)
    except Exception as e:
        print(f"[!] Could not read {file_path}: {e}")
        skipped.append(taxon)
        continue

    # Try to locate KO column by checking values
    ko_column = None
    for col in df.columns:
        sample_values = df[col].dropna().astype(str).head(10)
        if sample_values.str.contains("ko:K").any():
            ko_column = col
            break

    if not ko_column:
        print(f"[!] No KO column in {os.path.basename(file_path)}")
        skipped.append(taxon)
        continue

    kos = df[ko_column].dropna().astype(str)
    for ko_list in kos:
        for ko in ko_list.split(","):
            if ko.startswith("ko:K"):
                records.append({'Taxon': taxon, 'KO': ko.replace("ko:", "")})

if not records:
    print("❌ No KO terms were extracted from any file. Exiting.")
    exit()

# Save KO table
ko_df = pd.DataFrame(records).drop_duplicates()
ko_df.to_csv(os.path.join(output_folder, "taxon_to_ko.csv"), index=False)
print("📄 Saved KO mappings to: taxon_to_ko.csv")

# Save skipped taxa list
pd.DataFrame(skipped, columns=["Skipped_Taxon"]).to_csv(os.path.join(output_folder, "skipped_taxa.csv"), index=False)
print(f"📄 Saved list of skipped taxa ({len(skipped)}): skipped_taxa.csv")

# === 3. Map KOs to KEGG pathways ===
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

# Map KO → pathways
ko_to_pathways = {}
print("🔗 Mapping KO terms to KEGG pathways...")
for ko in sorted(ko_df["KO"].unique()):
    ko_to_pathways[ko] = get_pathways_for_ko(ko)
    time.sleep(0.5)

# Get pathway names
pathway_names = {}
print("🧠 Fetching pathway names...")
for pathway_list in ko_to_pathways.values():
    for pid in pathway_list:
        if pid not in pathway_names:
            pathway_names[pid] = get_pathway_name(pid)
            time.sleep(0.3)

# === 4. Final Output ===
final = []
for row in records:
    taxon, ko = row['Taxon'], row['KO']
    for pid in ko_to_pathways.get(ko, []):
        final.append({
            'Taxon': taxon,
            'KO': ko,
            'Pathway_ID': pid,
            'Pathway_Name': pathway_names.get(pid, "Unknown")
        })

final_df = pd.DataFrame(final).drop_duplicates()
final_df.to_csv(os.path.join(output_folder, "taxon_kegg_pathways.csv"), index=False)
print("🎉 Done! Pathway output saved to: taxon_kegg_pathways.csv")
