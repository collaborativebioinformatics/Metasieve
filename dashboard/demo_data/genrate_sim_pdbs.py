import os
import requests

output_dir = "metagenomic_dark_matter_pdbs"
os.makedirs(output_dir, exist_ok=True)

print("Fetching 500 unclassified / metagenomic dark matter IDs from UniProt...")

# Query UniProt for unreviewed proteins associated with metagenomes
uniprot_url = "https://rest.uniprot.org/uniprotkb/search?query=metagenome+AND+reviewed:false&size=500&fields=accession"
# query for strictly dark matter : uniprot_url = "https://rest.uniprot.org/uniprotkb/search?query=metagenome+AND+reviewed:false+AND+protein_name:%22uncharacterized+protein%22&size=500&fields=accession"
try:
    response = requests.get(uniprot_url, timeout=10)
    if response.status_code == 200:
        data = response.json()
        known_ids = [entry["primaryAccession"] for entry in data.get("results", [])]
        print(f"Successfully retrieved {len(known_ids)} dark matter IDs. Starting downloads...\n")
    else:
        print(f"Failed to fetch from UniProt API: Status {response.status_code}")
        known_ids = []
except Exception as e:
    print(f"Connection error while fetching IDs: {e}")
    known_ids = []

session = requests.Session()
success_count = 0
target_count = len(known_ids)

for uid in known_ids:
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uid}"
    try:
        # (Connect timeout, Read timeout) prevents infinite hanging
        response = session.get(url, timeout=(5, 10))
        
        if response.status_code == 200:
            data = response.json()
            if data and isinstance(data, list) and len(data) > 0:
                pdb_url = data[0].get("pdbUrl")
                if pdb_url:
                    pdb_response = session.get(pdb_url, timeout=(5, 10))
                    if pdb_response.status_code == 200:
                        out_path = os.path.join(output_dir, f"{uid}.pdb")
                        with open(out_path, "w", encoding="utf-8") as f:
                            f.write(pdb_response.text)
                        
                        success_count += 1
                        print(f"[{success_count}/{target_count}] Downloaded: {uid}.pdb")
                        
    except requests.exceptions.RequestException:
        # Silently skip network hiccups and move forward
        continue

print(f"\nDone! Successfully downloaded {success_count} metagenomic dark matter PDB files to '{output_dir}/'.")