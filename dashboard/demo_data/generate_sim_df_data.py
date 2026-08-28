import os
import glob
import random
import numpy as np

def parse_pdb_file(filepath):
    """Extracts basic metrics from a PDB file to ground the simulation."""
    atom_count = 0
    residues = set()
    plddt_scores = []
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_count += 1
                res_id = line[22:26].strip()
                residues.add(res_id)
                
                try:
                    b_factor = float(line[60:66].strip())
                    plddt_scores.append(b_factor)
                except ValueError:
                    pass
                    
    avg_plddt = np.mean(plddt_scores) if plddt_scores else float('nan')
    return {
        "filename": os.path.basename(filepath),
        "total_atoms": atom_count,
        "residue_count": len(residues),
        "mean_plddt": round(float(avg_plddt), 2)
    }

def generate_pipeline_simulations(pdb_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    
    total_files = len(pdb_files)
    if total_files == 0:
        print(f"No PDB files found in '{pdb_dir}'! Check that the folder contains .pdb files.")
        return

    print(f"Found {total_files} PDB files. Processing simulation reports with heterogeneous distribution...\n")

    # Shuffle files to randomly assign classified vs unclassified and pLDDT variations
    random.shuffle(pdb_files)
    split_index = int(total_files * 0.40)
    classified_files = set(pdb_files[:split_index]) # 40% classified

    # --- 1. KRAKEN 2 READ-LEVEL OUTPUT (.out format) ---
    kraken_out_path = os.path.join(output_dir, "simulated_kraken2_output.out")
    with open(kraken_out_path, "w") as k_out:
        for pdb in pdb_files:
            metrics = parse_pdb_file(pdb)
            seq_id = metrics["filename"].replace(".pdb", "")
            seq_length = max(50, metrics["residue_count"] * 3)
            
            if pdb in classified_files:
                # Classified entry (40%)
                tax_id = random.randint(1000, 99999)
                k_out.write(f"C\t{seq_id}\t{tax_id}\t{seq_length}\t12:15|15:10\n")
            else:
                # Unclassified dark matter entry (60%)
                k_out.write(f"U\t{seq_id}\t0\t{seq_length}\t0:12|1:15\n")

    # --- 2. SIMULATE SEQSCREEN RESULTS ---
    seqscreen_path = os.path.join(output_dir, "simulated_seqscreen_results.txt")
    with open(seqscreen_path, "w") as s_out:
        s_out.write("sequence_id\tfunctional_category\tec_number\tpathogenicity_score\tthreat_level\tnovelty_status\n")
        
        for pdb in pdb_files:
            metrics = parse_pdb_file(pdb)
            seq_id = metrics["filename"].replace(".pdb", "")
            
            if pdb in classified_files:
                cat = "Enzyme Metabolism"
                ec = f"{random.randint(1,6)}.{random.randint(1,5)}.{random.randint(1,9)}.{random.randint(1,15)}"
                score = round(random.uniform(0.0, 0.05), 3)
                threat = "Low Risk"
                novelty = "Known / Annotated"
            else:
                cat = random.choices(
                    ["Uncharacterized Dark Matter", "Hypothetical Protein"], 
                    weights=[0.8, 0.2]
                )[0]
                ec = "-"
                score = round(random.uniform(0.0, 0.2), 3)
                threat = "Indeterminate/Dark Matter"
                novelty = "Novel/Dark Matter"
            
            s_out.write(f"{seq_id}\t{cat}\t{ec}\t{score}\t{threat}\t{novelty}\n")

    # --- 3. SIMULATE ESMFOLD / STRUCTURAL METRICS LOG ---
    esm_path = os.path.join(output_dir, "simulated_esmfold_summary.csv")
    with open(esm_path, "w") as e_out:
        e_out.write("sequence_id,total_atoms,residue_count,mean_plddt,pTM,inferred_status\n")
        
        for i, pdb in enumerate(pdb_files):
            metrics = parse_pdb_file(pdb)
            seq_id = metrics["filename"].replace(".pdb", "")
            
            # Force a balanced variety of pLDDT < 70 and > 70 across the dataset
            # Alternating or randomizing to ensure robust graphing options
            if i % 2 == 0:
                forced_plddt = round(random.uniform(40.0, 69.5), 2)
                status = "Low Confidence / Structural Dark Matter"
            else:
                forced_plddt = round(random.uniform(70.5, 95.0), 2)
                status = "Moderate/High Confidence"
                
            ptm = round(forced_plddt / 100.0 * random.uniform(0.8, 1.0), 3)
            
            e_out.write(f"{seq_id},{metrics['total_atoms']},{metrics['residue_count']},{forced_plddt},{ptm},{status}\n")

    print(f"Heterogeneous simulation complete! Processed {total_files} files.")
    print(f" - Kraken Output (40% Classified / 60% Unclassified): {kraken_out_path}")
    print(f" - SeqScreen Table: {seqscreen_path}")
    print(f" - ESMFold Metrics (Mixed pLDDT <70 and >70): {esm_path}")

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    target_pdb_dir = os.path.join(script_dir, "metagenomic_dark_matter_pdbs")
    target_out_dir = os.path.join(script_dir, "pipeline_outputs")
    
    generate_pipeline_simulations(pdb_dir=target_pdb_dir, output_dir=target_out_dir)