#!/usr/bin/env python3
"""
Script to fetch human cytoplasmic ribosomal protein information from UniProt.
Creates a comprehensive table of all human ribosomal proteins with UniProt IDs and unified nomenclature.
"""

import requests
import time
import csv
import json

# Complete list of human ribosomal protein genes based on Wikipedia and literature
# Small subunit (40S) proteins
SMALL_SUBUNIT = {
    'RPSA': 'uS2',
    'RPS3': 'uS3',
    'RPS9': 'uS4',
    'RPS2': 'uS5',
    'RPS6': 'eS6',
    'RPS5': 'uS7',
    'RPS8': 'eS8',
    'RPS16': 'uS9',
    'RPS20': 'uS10',
    'RPS14': 'uS11',
    'RPS23': 'uS12',
    'RPS18': 'uS13',
    'RPS29': 'uS14',
    'RPS13': 'uS15',
    'RPS11': 'uS17',
    'RPS15': 'uS19',
    'RPS19': 'eS19',
    'RPS21': 'eS21',
    'RPS24': 'eS24',
    'RPS25': 'eS25',
    'RPS26': 'eS26',
    'RPS27': 'eS27',
    'RPS28': 'eS28',
    'RPS30': 'eS30',
    'RPS27A': 'eS31',
    'RACK1': 'RACK1',
    # Additional RPS genes
    'RPS3A': 'eS1',
    'RPS7': 'eS7',
    'RPS10': 'eS10',
    'RPS12': 'eS12',
    'RPS15A': 'uS8',
    'RPS17': 'eS17',
    'RPS4X': 'eS4',
    'RPS4Y1': 'eS4',
    'RPS4Y2': 'eS4',
}

# Large subunit (60S) proteins
LARGE_SUBUNIT = {
    'RPL10A': 'uL1',
    'RPL8': 'uL2',
    'RPL3': 'uL3',
    'RPL4': 'uL4',
    'RPL11': 'uL5',
    'RPL6': 'eL6',
    'RPL7A': 'eL8',
    'RPLP0': 'uL10',
    'RPL12': 'uL11',
    'RPL13A': 'uL13',
    'RPL13': 'eL13',
    'RPL23': 'uL14',
    'RPL14': 'eL14',
    'RPL27A': 'uL15',
    'RPL15': 'eL15',
    'RPL10': 'uL16',
    'RPL5': 'uL18',
    'RPL18': 'eL18',
    'RPL19': 'eL19',
    'RPL18A': 'eL20',
    'RPL21': 'eL21',
    'RPL17': 'uL22',
    'RPL22': 'eL22',
    'RPL23A': 'uL23',
    'RPL26': 'uL24',
    'RPL24': 'eL24',
    'RPL27': 'eL27',
    'RPL28': 'eL28',
    'RPL35': 'uL29',
    'RPL29': 'eL29',
    'RPL7': 'uL30',
    'RPL30': 'eL30',
    'RPL31': 'eL31',
    'RPL32': 'eL32',
    'RPL35A': 'eL33',
    'RPL34': 'eL34',
    'RPL36': 'eL36',
    'RPL37': 'eL37',
    'RPL38': 'eL38',
    'RPL39': 'eL39',
    'RPL40': 'eL40',
    'RPL41': 'eL41',
    'RPL36A': 'eL42',
    'RPL37A': 'eL43',
    'RPLP1': 'P1',
    'RPLP2': 'P2',
    # Additional RPL genes
    'RPL9': 'uL6',
    'RPL22L1': 'eL22-like',
    'RPL26L1': 'uL24-like',
}

# Special cases
SPECIAL_CASES = {
    'FAU': 'eS30-ubiquitin',  # ubiquitin-RPL40 precursor
}

def query_uniprot(gene_name, organism='9606'):
    """
    Query UniProt for a specific gene in human (taxonomy 9606).
    Returns the primary UniProt accession and protein name.
    """
    base_url = 'https://rest.uniprot.org/uniprotkb/search'

    # Query for the gene in humans, reviewed (Swiss-Prot) entries only
    query = f'(gene:{gene_name}) AND (organism_id:{organism}) AND (reviewed:true)'

    params = {
        'query': query,
        'format': 'json',
        'size': 5  # Get top 5 results to check
    }

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()

        if data.get('results'):
            # Get the first (primary) entry
            entry = data['results'][0]
            accession = entry['primaryAccession']
            protein_name = entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '')

            # Also get gene names to verify
            gene_names = entry.get('genes', [])
            primary_gene = gene_names[0]['geneName']['value'] if gene_names else ''

            return {
                'accession': accession,
                'protein_name': protein_name,
                'primary_gene': primary_gene
            }
        else:
            return None

    except Exception as e:
        print(f"Error querying {gene_name}: {e}")
        return None

def fetch_all_ribosomal_proteins():
    """
    Fetch UniProt information for all human cytoplasmic ribosomal proteins.
    """
    all_proteins = []

    print("Fetching Small Subunit (40S) proteins...")
    for gene, unified_name in sorted(SMALL_SUBUNIT.items()):
        print(f"  Querying {gene}...")
        result = query_uniprot(gene)
        if result:
            all_proteins.append({
                'Gene_Symbol': gene,
                'Unified_Nomenclature': unified_name,
                'UniProt_ID': result['accession'],
                'Protein_Name': result['protein_name'],
                'Subunit': '40S',
                'Verified_Gene': result['primary_gene']
            })
        else:
            print(f"    WARNING: No result found for {gene}")
            all_proteins.append({
                'Gene_Symbol': gene,
                'Unified_Nomenclature': unified_name,
                'UniProt_ID': 'NOT_FOUND',
                'Protein_Name': '',
                'Subunit': '40S',
                'Verified_Gene': ''
            })
        time.sleep(0.2)  # Be polite to UniProt API

    print("\nFetching Large Subunit (60S) proteins...")
    for gene, unified_name in sorted(LARGE_SUBUNIT.items()):
        print(f"  Querying {gene}...")
        result = query_uniprot(gene)
        if result:
            all_proteins.append({
                'Gene_Symbol': gene,
                'Unified_Nomenclature': unified_name,
                'UniProt_ID': result['accession'],
                'Protein_Name': result['protein_name'],
                'Subunit': '60S',
                'Verified_Gene': result['primary_gene']
            })
        else:
            print(f"    WARNING: No result found for {gene}")
            all_proteins.append({
                'Gene_Symbol': gene,
                'Unified_Nomenclature': unified_name,
                'UniProt_ID': 'NOT_FOUND',
                'Protein_Name': '',
                'Subunit': '60S',
                'Verified_Gene': ''
            })
        time.sleep(0.2)

    print("\nFetching special cases...")
    for gene, unified_name in sorted(SPECIAL_CASES.items()):
        print(f"  Querying {gene}...")
        result = query_uniprot(gene)
        if result:
            all_proteins.append({
                'Gene_Symbol': gene,
                'Unified_Nomenclature': unified_name,
                'UniProt_ID': result['accession'],
                'Protein_Name': result['protein_name'],
                'Subunit': 'Special',
                'Verified_Gene': result['primary_gene']
            })
        else:
            print(f"    WARNING: No result found for {gene}")
            all_proteins.append({
                'Gene_Symbol': gene,
                'Unified_Nomenclature': unified_name,
                'UniProt_ID': 'NOT_FOUND',
                'Protein_Name': '',
                'Subunit': 'Special',
                'Verified_Gene': ''
            })
        time.sleep(0.2)

    return all_proteins

def save_to_csv(proteins, filename='human_ribosomal_proteins.csv'):
    """
    Save the protein data to a CSV file.
    """
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        fieldnames = ['Gene_Symbol', 'Unified_Nomenclature', 'UniProt_ID',
                     'Protein_Name', 'Subunit', 'Verified_Gene']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(proteins)

    print(f"\nData saved to {filename}")

def print_summary(proteins):
    """
    Print a summary of the fetched data.
    """
    total = len(proteins)
    found = sum(1 for p in proteins if p['UniProt_ID'] != 'NOT_FOUND')
    small = sum(1 for p in proteins if p['Subunit'] == '40S')
    large = sum(1 for p in proteins if p['Subunit'] == '60S')

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total proteins queried: {total}")
    print(f"Successfully found: {found}")
    print(f"Not found: {total - found}")
    print(f"Small subunit (40S): {small}")
    print(f"Large subunit (60S): {large}")
    print("="*60)

if __name__ == '__main__':
    print("="*60)
    print("Human Cytoplasmic Ribosomal Protein UniProt Fetcher")
    print("="*60)

    proteins = fetch_all_ribosomal_proteins()
    save_to_csv(proteins)
    print_summary(proteins)

    print("\nDone!")
