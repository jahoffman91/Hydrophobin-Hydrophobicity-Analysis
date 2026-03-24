import pandas as pd
import plotly.graph_objects as go
from Bio import SeqIO
import argparse
import os

# Define genetic code dictionary at module level
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def read_fasta(fasta_file):
    """
    Read sequences from a FASTA file
    Returns a list of tuples (sequence_id, sequence)
    """
    sequences = []
    try:
        # Get absolute path and check if file exists
        abs_path = os.path.abspath(fasta_file)
        print(f"Looking for file at: {abs_path}")
        print(f"File exists: {os.path.exists(abs_path)}")
        print(f"Current working directory: {os.getcwd()}")
        print(f"Directory contents: {os.listdir('.')}")
        
        with open(fasta_file, 'r') as f:
            print("File contents:")
            content = f.read()
            print(content)
            if not content.strip():
                print("Warning: File is empty!")
            f.seek(0)
            
            for record in SeqIO.parse(f, "fasta"):
                print(f"Found sequence: {record.id}")
                sequences.append((record.id, str(record.seq)))
                
        if not sequences:
            print("No sequences found. Check if file is in correct FASTA format:")
            print("Expected format:")
            print(">sequence_name")
            print("SEQUENCE")
            
    except FileNotFoundError:
        print(f"Error: Could not find file {fasta_file}")
    except Exception as e:
        print(f"Error reading file: {str(e)}")
    
    return sequences

def analyze_sequence_hydrophobicity(sequence, hydrophobicity_scales):
    """
    Analyze a sequence using all hydrophobicity scales
    Returns a DataFrame with hydrophobicity scores for each position
    """
    # Initialize results dictionary
    results = {scale: [] for scale in hydrophobicity_scales.columns}
    
    # Get hydrophobicity score for each position
    for aa in sequence:
        for scale in hydrophobicity_scales.columns:
            try:
                # Convert to float to ensure numeric value
                score = float(hydrophobicity_scales.loc[aa, scale])
                results[scale].append(score)
            except (KeyError, ValueError) as e:
                print(f"Warning: Issue with amino acid {aa} in scale {scale}: {e}")
                results[scale].append(None)
    
    # Create DataFrame with results and ensure numeric type
    df_results = pd.DataFrame(results).astype(float)
    df_results.index = [f"{i+1}:{aa}" for i, aa in enumerate(sequence)]
    return df_results

def translate_dna_and_track_uppercase(dna_sequence):
    """
    Translates DNA sequence and tracks which amino acids came from uppercase codons
    Returns: (amino acid sequence, list of positions with uppercase codons)
    """
    amino_acids = []
    uppercase_positions = []
    
    # Process sequence in codons
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:  # Only process complete codons
            is_uppercase = codon.isupper()
            codon_upper = codon.upper()
            
            if codon_upper in GENETIC_CODE:
                aa = GENETIC_CODE[codon_upper]
                if aa != '*':  # Skip stop codons
                    amino_acids.append(aa)
                    if is_uppercase:
                        uppercase_positions.append(len(amino_acids) - 1)
    
    return ''.join(amino_acids), uppercase_positions

def plot_hydrophobicity_profile(df_results, title="Hydrophobicity Profile", aa_per_row=50, 
                              uppercase_positions=None, signal_sequence_positions=None):
    # Calculate mean and std of hydrophobicity for each position
    mean_hydrophobicity = df_results.mean(axis=1)
    std_hydrophobicity = df_results.std(axis=1)
    max_hydrophobicity = mean_hydrophobicity.max()  # Get the maximum hydrophobicity value
    
    # Calculate 3-residue rolling average
    positions = list(range(1, len(mean_hydrophobicity) + 1))
    amino_acids = [idx.split(':')[1] for idx in df_results.index]
    rolling_avg = pd.Series(mean_hydrophobicity).rolling(window=3, center=True).mean()
    
    # Calculate number of rows needed
    n_rows = (len(positions) - 1) // aa_per_row + 1
    
    # Create subplots
    fig = go.Figure()
    
    # Add signal sequence highlight before the row loop
    if signal_sequence_positions and signal_sequence_positions:
        # Get the start and end of signal sequence
        start_pos = signal_sequence_positions[0] + 1
        end_pos = signal_sequence_positions[-1] + 1
        
        fig.add_shape(
            type="rect",
            x0=start_pos-0.5,
            x1=end_pos+0.5,
            y0=0,
            y1=100,
            fillcolor="rgba(0,0,255,0.1)",
            layer="below",
            line_width=0,
            yref="y",
            xref="x"
        )
        
        # Add single annotation for the signal sequence
        fig.add_annotation(
            x=start_pos,
            y=100,
            text="Signal Sequence",
            showarrow=False,
            yshift=10,
            xref="x",
            yref="y"
        )
    
    for row in range(n_rows):
        start_idx = row * aa_per_row
        end_idx = min((row + 1) * aa_per_row, len(positions))
        
        # Get data for this row
        row_positions = positions[start_idx:end_idx]
        row_aa = amino_acids[start_idx:end_idx]
        row_mean = mean_hydrophobicity[start_idx:end_idx]
        row_rolling = rolling_avg[start_idx:end_idx]
        
        # Add highlight bars for uppercase positions (red background)
        if uppercase_positions:
            for pos in row_positions:
                if (pos-1) in uppercase_positions:
                    fig.add_vrect(
                        x0=pos-0.5,
                        x1=pos+0.5,
                        fillcolor="rgba(255,0,0,0.1)",  # Light red
                        layer="below",
                        line_width=0,
                        xref=f'x{row+1}' if row > 0 else 'x',
                        yref=f'y{row+1}' if row > 0 else 'y',
                        showlegend=False  # Don't show in legend
                    )
        
        hover_text = []
        for pos, aa, mean_h in zip(row_positions, row_aa, row_mean):
            text = (
                f"Position: {pos}<br>"
                f"Amino Acid: {aa}<br>"
            )
            if signal_sequence_positions:
                is_signal = "YES" if (pos-1) in signal_sequence_positions else "NO"
                text += f"Signal Sequence: {is_signal}<br>"
            if uppercase_positions:
                is_uppercase = "YES" if (pos-1) in uppercase_positions else "NO"
                text += f"Uppercase Codon: {is_uppercase}<br>"
            text += f"Hydrophobicity: {mean_h:.2f}"
            hover_text.append(text)
        
        # Add mean hydrophobicity points with color based on uppercase status
        marker_colors = ['red' if (pos-1) in (uppercase_positions or []) else 'blue' 
                        for pos in row_positions]
        
        # Add max hydrophobicity line
        fig.add_trace(go.Scatter(
            x=row_positions,
            y=[max_hydrophobicity] * len(row_positions),
            mode='lines',
            name='Max Hydrophobicity',
            line=dict(color='gray', width=1, dash='dot'),
            xaxis=f'x{row+1}' if row > 0 else 'x',
            yaxis=f'y{row+1}' if row > 0 else 'y',
            showlegend=row == 0  # Only show in legend for first row
        ))
        
        # Add points
        fig.add_trace(go.Scatter(
            x=row_positions,
            y=row_mean,
            mode='markers',
            name='Hydrophobicity Score',
            marker=dict(color=marker_colors, size=6),
            hovertext=hover_text,
            hoverinfo='text',
            xaxis=f'x{row+1}' if row > 0 else 'x',
            yaxis=f'y{row+1}' if row > 0 else 'y',
            showlegend=row == 0
        ))
        
        # Add rolling average line
        fig.add_trace(go.Scatter(
            x=row_positions,
            y=row_rolling,
            mode='lines',
            name='3-residue average',
            line=dict(color='red', width=2),
            xaxis=f'x{row+1}' if row > 0 else 'x',
            yaxis=f'y{row+1}' if row > 0 else 'y',
            showlegend=row == 0
        ))
        
        # Create x-axis labels for this row
        x_labels = []
        for pos, aa in zip(row_positions, row_aa):
            if (pos-1) in (uppercase_positions or []):
                label = f"{pos}<br><span style='color: red; font-weight: bold'>{aa}</span>"
            else:
                label = f"{pos}<br>{aa}"
            x_labels.append(label)
        
        # Update axes for this row
        xaxis_key = f'xaxis{row+1}' if row > 0 else 'xaxis'
        yaxis_key = f'yaxis{row+1}' if row > 0 else 'yaxis'
        
        # Calculate evenly spaced domains with gaps between rows
        spacing = 0.1
        row_height = (1.0 - (n_rows - 1) * spacing) / n_rows
        y_domain = [
            max(0, 1 - (row + 1) * row_height - row * spacing),
            min(1, 1 - row * row_height - row * spacing)
        ]
        
        # Update layout for this row
        fig.update_layout(**{
            xaxis_key: dict(
                ticktext=x_labels,
                tickvals=row_positions,
                tickmode='array',
                tickangle=0,
                tickfont=dict(size=8),
                title='Sequence Position' if row == n_rows-1 else '',
                range=[min(row_positions)-0.5, 
                      max(row_positions)+0.5 if row != n_rows-1 
                      else start_idx + aa_per_row + 0.5]
            ),
            yaxis_key: dict(
                range=[0, 100],
                title='Hydrophobicity Score (0-100)' if row == 0 else '',
                zeroline=False,
                domain=y_domain
            )
        })
    
    # Update overall layout
    fig.update_layout(
        title=title,
        height=300 * n_rows,
        showlegend=True,
        grid=dict(rows=n_rows, columns=1, pattern='independent'),
    )
    
    return fig

def get_uppercase_positions_from_highlighted(highlighted_fasta, sequences):
    """
    Get the positions of amino acids that correspond to uppercase codons,
    excluding signal sequences
    Returns: (uppercase_map, trimmed_sequences)
    """
    uppercase_map = {}
    trimmed_sequences = {}  # Store sequences with signal sequences removed
    highlighted_seqs = read_fasta(highlighted_fasta)
    
    # Create a map of full sequences for alignment
    full_seq_map = {seq_id: seq for seq_id, seq in sequences}
    
    for seq_id, dna_seq in highlighted_seqs:
        print(f"\nProcessing highlighted sequence {seq_id}")
        
        if seq_id not in full_seq_map:
            print(f"Warning: Sequence {seq_id} not found in full sequences")
            continue
            
        # Get the full protein sequence
        full_seq = full_seq_map[seq_id]
        
        # Process sequence and get translated sequence
        translated_seq = ''
        uppercase_positions = []
        
        # Process sequence in codons and track uppercase positions
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i+3]
            if len(codon) == 3:
                is_uppercase = codon.isupper()
                codon_upper = codon.upper()
                if codon_upper in GENETIC_CODE:
                    aa = GENETIC_CODE[codon_upper]
                    if aa != '*':
                        translated_seq += aa
                        if is_uppercase:
                            pos = len(translated_seq) - 1
                            uppercase_positions.append(pos)
        
        # Find where translated sequence aligns in full sequence
        if translated_seq in full_seq:
            offset = full_seq.index(translated_seq)
            # Store the mature sequence (without signal sequence)
            trimmed_sequences[seq_id] = full_seq[offset:]
            # Adjust positions (no need to add offset since we're trimming the sequence)
            uppercase_map[seq_id] = uppercase_positions
            print(f"Trimmed signal sequence of length {offset} from {seq_id}")
        else:
            print(f"Warning: Could not align highlighted sequence for {seq_id}")
    
    return uppercase_map, trimmed_sequences

def analyze_sequences(fasta_file, highlighted_fasta, output_prefix="hydrophobicity_analysis"):
    """
    Analyze sequences from a FASTA file and generate plots
    """
    print(f"Loading scales from hydrophobicity_scales_normalized.csv...")
    scales_df = pd.read_csv('hydrophobicity_scales_normalized.csv', index_col=0)
    hydrophobicity_scales = scales_df.iloc[3:].astype(float)
    print(f"Loaded {len(hydrophobicity_scales.columns)} scales")
    
    print(f"Reading sequences from {fasta_file}...")
    sequences = read_fasta(fasta_file)
    print(f"Found {len(sequences)} sequences")
    
    # Get uppercase positions and trimmed sequences
    uppercase_map, trimmed_sequences = get_uppercase_positions_from_highlighted(highlighted_fasta, sequences)
    
    # Process each sequence
    for seq_id, seq in trimmed_sequences.items():
        print(f"\nProcessing sequence: {seq_id}")
        
        # Get uppercase positions for this sequence
        uppercase_pos = uppercase_map.get(seq_id, [])
        
        # Analyze sequence
        results = analyze_sequence_hydrophobicity(seq, hydrophobicity_scales)
        
        # Create plot with highlighting
        print(f"Creating plot...")
        fig = plot_hydrophobicity_profile(
            results, 
            f"{seq_id} Hydrophobicity Profile",
            uppercase_positions=uppercase_pos
        )
        
        # Save results
        print(f"Saving results to {output_prefix}_{seq_id}.csv")
        results.to_csv(f"{output_prefix}_{seq_id}.csv")
        print(f"Saving plot to {output_prefix}_{seq_id}.html")
        fig.write_html(f"{output_prefix}_{seq_id}.html")
        
        # Print summary statistics
        print(f"\nSequence: {seq_id}")
        print(f"Length: {len(seq)}")
        print("\nHydrophobicity Statistics:")
        print(results.describe())
        print("\n" + "="*50)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze protein sequence hydrophobicity from FASTA file')
    parser.add_argument('fasta_file', help='Input FASTA file containing protein sequences')
    parser.add_argument('highlighted_fasta', help='FASTA file with uppercase/lowercase highlighting')
    parser.add_argument('--output', '-o', default='hydrophobicity_analysis',
                      help='Prefix for output files (default: hydrophobicity_analysis)')
    
    args = parser.parse_args()
    
    analyze_sequences(args.fasta_file, args.highlighted_fasta, args.output)
