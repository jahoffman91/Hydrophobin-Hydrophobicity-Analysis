import pandas as pd
import numpy as np
import plotly.graph_objects as go

def determine_scale_direction(df):
    """
    Determine if a scale is positive-hydrophobic or negative-hydrophobic
    by comparing average values of known hydrophobic vs hydrophilic amino acids
    """
    # Define known hydrophobic and hydrophilic amino acids
    hydrophobic_aas = ['I', 'L', 'W', 'V', 'F']  # Most hydrophobic
    hydrophilic_aas = ['R', 'K', 'D', 'E']       # Most hydrophilic
    
    scale_directions = {}
    
    for column in df.columns:
        # Calculate average value for hydrophobic and hydrophilic AAs
        hydrophobic_avg = df.loc[hydrophobic_aas, column].mean()
        hydrophilic_avg = df.loc[hydrophilic_aas, column].mean()
        
        # If hydrophobic AAs have higher values, scale is positive-hydrophobic
        # If hydrophobic AAs have lower values, scale is negative-hydrophobic
        scale_directions[column] = 1 if hydrophobic_avg > hydrophilic_avg else -1
        
        print(f"\n{column}:")
        print(f"Hydrophobic average: {hydrophobic_avg:.2f}")
        print(f"Hydrophilic average: {hydrophilic_avg:.2f}")
        print(f"Direction: {'positive' if scale_directions[column] > 0 else 'negative'}-hydrophobic")
    
    return scale_directions

def normalize_scales(df, scale_directions):
    """
    Normalize each scale to 0-100 range where 100 is most hydrophobic
    """
    df_normalized = pd.DataFrame(index=df.index)
    
    for column in df.columns:
        values = df[column].values
        direction = scale_directions[column]
        
        if direction == -1:
            # If negative-hydrophobic, invert the values
            values = -values
            
        # Normalize to 0-100 range
        min_val = np.min(values)
        max_val = np.max(values)
        normalized = (values - min_val) * 100 / (max_val - min_val)
        
        df_normalized[column] = normalized
        
    return df_normalized

def create_hydrophobicity_boxplot(df_normalized):
    """
    Create a box plot of hydrophobicity scores for each amino acid
    """
    # Get author names from the first row
    authors = df_normalized.iloc[0]
    # Create a dictionary mapping column names to author names
    author_dict = {col: (author if not pd.isna(author) else col) 
                  for col, author in zip(authors.index, authors)}
    
    # Get the data, skipping metadata rows
    data = df_normalized.iloc[3:].astype(float)
    
    # AA code mapping
    aa_codes = {
        'A': 'Ala (A)', 'C': 'Cys (C)', 'D': 'Asp (D)', 'E': 'Glu (E)',
        'F': 'Phe (F)', 'G': 'Gly (G)', 'H': 'His (H)', 'I': 'Ile (I)',
        'K': 'Lys (K)', 'L': 'Leu (L)', 'M': 'Met (M)', 'N': 'Asn (N)',
        'P': 'Pro (P)', 'Q': 'Gln (Q)', 'R': 'Arg (R)', 'S': 'Ser (S)',
        'T': 'Thr (T)', 'V': 'Val (V)', 'W': 'Trp (W)', 'Y': 'Tyr (Y)'
    }
    
    # Calculate median for each amino acid and sort in descending order
    medians = data.median(axis=1)
    sorted_aas = medians.sort_values(ascending=False).index
    
    # Create box plot
    fig = go.Figure()
    
    # Default plotly colors
    colors = [
        '#1f77b4',  # blue
        '#ff7f0e',  # orange
        '#2ca02c',  # green
        '#d62728',  # red
        '#9467bd',  # purple
        '#8c564b',  # brown
        '#e377c2',  # pink
        '#7f7f7f',  # gray
        '#bcbd22',  # yellow-green
        '#17becf'   # cyan
    ]
    
    # For each amino acid (in sorted order), create a box plot
    for i, aa in enumerate(sorted_aas):
        # Create hover text with author names and values
        hover_text = []
        for col, value in data.loc[aa].items():
            author = author_dict[col]  # Use the dictionary instead of list indexing
            hover_text.append(f"Author: {author}<br>Score: {value:.2f}")
        
        fig.add_trace(go.Box(
            y=data.loc[aa],
            name=aa_codes[aa],
            boxpoints='all',
            jitter=0.3,
            pointpos=0,
            marker=dict(
                size=4,
                color='black'
            ),
            boxmean=True,
            line=dict(
                color='black',
                width=2
            ),
            fillcolor=colors[i % len(colors)],
            showlegend=False,
            hovertext=hover_text,
            hoverinfo='text'
        ))
    
    # Add invisible scatter traces just for the legend
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='black', width=2),
        name='Median'
    ))
    
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='black', width=2, dash='dash'),
        name='Mean'
    ))
    
    # Update layout
    fig.update_layout(
        title='Hydrophobicity Scores Distribution by Amino Acid',
        yaxis_title='Normalized Hydrophobicity Score (0-100)',
        xaxis_title='Amino Acid',
        height=600,
        width=1000,
        yaxis=dict(range=[0, 110]),
        showlegend=True
    )
    
    return fig

if __name__ == "__main__":
    # Read the CSV file
    df = pd.read_csv('hydrophobicity_scales.csv', index_col=0)
    
    # Extract just the amino acid rows (skip metadata)
    df_values = df.iloc[3:]
    df_values = df_values.astype(float)
    
    # Determine scale directions
    scale_directions = determine_scale_direction(df_values)
    
    # Normalize scales
    df_normalized = normalize_scales(df_values, scale_directions)
    
    print("\nNormalized values (0-100 scale):")
    print(df_normalized)
    
    # Save normalized values
    output_file = 'hydrophobicity_scales_normalized.csv'
    
    # Combine metadata with normalized values
    df_output = pd.concat([df.iloc[:3], df_normalized])
    df_output.to_csv(output_file)
    print(f"\nNormalized scales saved to {output_file}")
    
    # After creating df_output
    print("Creating box plot visualization...")
    fig = create_hydrophobicity_boxplot(df_output)
    fig.write_html('hydrophobicity_distribution.html')
    print("Box plot saved as hydrophobicity_distribution.html")