import requests
from bs4 import BeautifulSoup
import pandas as pd
import time

def get_scale_metadata(content):
    """Extract metadata from the page content"""
    metadata = {}
    
    # Extract title
    if "Title:" in content:
        title_section = content.split("Title:")[1].split("Author(s):")[0]
        metadata['Title'] = title_section.strip()
    
    # Extract authors
    if "Author(s):" in content:
        authors_section = content.split("Author(s):")[1].split("Reference:")[0]
        metadata['Authors'] = authors_section.strip()
    
    # Extract reference
    if "Reference:" in content:
        ref_section = content.split("Reference:")[1].split("Amino acid scale values:")[0]
        metadata['Reference'] = ref_section.strip()
    
    return metadata

def get_scale_data(url):
    """Get hydrophobicity values and metadata from a specific scale page"""
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')
    content = soup.get_text()
    
    # Get metadata
    metadata = get_scale_metadata(content)
    
    # Dictionary for three to one letter conversion
    aa_dict = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
        'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
        'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
        'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
        'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }
    
    # Get scale values
    scale_data = {}
    if "Amino acid scale values:" in content:
        values_section = content.split("Amino acid scale values:")[1].split("SIB logo")[0]
        
        for line in values_section.strip().split('\n'):
            line = line.strip()
            if line and ':' in line:
                aa, value = line.split(':')
                aa = aa.strip()
                if len(aa) == 3 and aa in aa_dict:  # Three letter code
                    try:
                        one_letter = aa_dict[aa]
                        value = float(value.strip())
                        scale_data[one_letter] = value
                    except ValueError:
                        continue
    
    return metadata, scale_data

def scrape_hydrophobicity_scales():
    base_url = "https://web.expasy.org/protscale/pscale"
    
    # Complete list of hydrophobicity scales from the website
    scales = [
        "Hphob.Eisenberg",      # Eisenberg et al.
        "Hphob.Sweet",          # Sweet et al. (OMH)
        "Hphob.Woods",          # Hopp & Woods
        "Hphob.Manavalan",      # Manavalan et al.
        "Hphob.Leo",            # Abraham & Leo
        "Hphob.Black",          # Black
        "Hphob.Breese",         # Bull & Breese
        "Hphob.Fauchere",       # Fauchere et al.
        "Hphob.Guy",            # Guy
        "Hphob.Janin",          # Janin
        "Hphob.Miyazawa",       # Miyazawa et al.
        "Hphob.Argos",          # Rao & Argos
        "Hphob.Roseman",        # Roseman
        "Hphob.Tanford",        # Tanford
        "Hphob.Wolfenden",      # Wolfenden et al.
        "Hphob.Welling",        # Welling & al
        "Hphob.Wilson",         # Wilson & al (HPLC)
        "Hphob.Parker",         # Parker & al (HPLC)
        "Hphob.pH3.4",          # Cowan (HPLC pH3.4)
        "Hphob.pH7.5",          # Cowan (HPLC pH7.5)
        "Hphob.mobility",       # Rf mobility
        "Hphob.Chothia",        # Chothia
        "Hphob.Rose"            # Rose & al
    ]
    
    # Dictionaries to store all data
    all_data = {}
    
    for scale in scales:
        print(f"Scraping {scale}...")
        scale_url = f"{base_url}/{scale}.html"
        metadata, scale_data = get_scale_data(scale_url)
        
        if scale_data:
            # Combine metadata and scale data
            combined_data = {
                'Title': metadata.get('Title', ''),
                'Authors': metadata.get('Authors', ''),
                'Reference': metadata.get('Reference', ''),
                **scale_data  # Add all amino acid values
            }
            all_data[scale] = combined_data
        
        time.sleep(1)  # Be nice to the server
    
    # Create DataFrame with metadata rows first, then amino acids
    df = pd.DataFrame(all_data)
    
    # Define the desired order of rows
    metadata_rows = ['Title', 'Authors', 'Reference']
    
    # Specified hydrophobic order for amino acids
    ordered_aas = [
        'I',  # Isoleucine
        'L',  # Leucine
        'W',  # Tryptophan
        'V',  # Valine
        'Y',  # Tyrosine
        'F',  # Phenylalanine
        'M',  # Methionine
        'C',  # Cysteine
        'A',  # Alanine
        'T',  # Threonine
        'H',  # Histidine
        'G',  # Glycine
        'S',  # Serine
        'Q',  # Glutamine
        'N',  # Asparagine
        'E',  # Glutamic acid
        'D',  # Aspartic acid
        'K',  # Lysine
        'R',  # Arginine
        'P'   # Proline
    ]
    
    desired_order = metadata_rows + ordered_aas
    df = df.reindex(desired_order)
    
    # Save to CSV
    output_file = 'hydrophobicity_scales.csv'
    df.to_csv(output_file)
    print(f"\nData saved to {output_file}")
    
    return df



if __name__ == "__main__":
    print("Starting hydrophobicity scale scraping from ExPASy ProtScale...")
    df = scrape_hydrophobicity_scales()
    print("\nAll rows of the data:")
    pd.set_option('display.max_rows', None)  # Show all rows
    print(df)

    # After creating the DataFrame:
    df_values = df.iloc[3:]  # Skip metadata rows
    df_values = df_values.astype(float)  # Convert to numeric
    
