"""
CellPyAbility_synergy.py is intended for dose combination analysis of 59 combinations and one vehicle.
This script should remain in the same directory as the other CellPyAbility scripts.
For more information, please see the README at https://github.com/bindralab/CellPyAbility.
"""

import tkinter as tk
from tkinter import ttk, filedialog

import toolbox as tb
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# Initialize toolbox
logger, base_dir = tb.logger, tb.base_dir

# Establish the GUI for experiment info
def synergy_gui():
    image_dir = ''
    image_dir_2 = ''

    def select_image_dir():
        nonlocal image_dir #opens a folder chooser
        image_dir = filedialog.askdirectory()

    def select_image_dir_2():
        nonlocal image_dir_2
        image_dir_2 = filedialog.askdirectory()

    # Create main GUI window
    root = tk.Tk()
    root.title('synergy')

    # Main Container Frame
    main_frame = ttk.Frame(root)
    main_frame.pack(padx = 20, pady = 20)

    # Left and right frames for side-by-side layout
    left_frame = ttk.Frame(main_frame)
    right_frame = ttk.Frame(main_frame)

    left_frame.grid(row=0, column=0, padx=15)
    right_frame.grid(row=0, column=1, padx=15)

    # Genotype 1
    ttk.Label(left_frame, text="Genotype 1 Parameters", font=('Arial', 12, 'bold')).pack()

    # Genotype 1 Fields
    entries = {}
    fields = [
        ('title_name', 'Enter the title of the experiment:'),
        ('x_drug', 'Enter the drug name for the horizontal gradient:'),
        ('x_top_conc', 'Enter the horizontal top concentration (uM):'),
        ('x_dilution', 'Enter the horizontal dilution factor (x-fold):'),
        ('y_drug', 'Enter the drug name for the vertical gradient:'),
        ('y_top_conc', 'Enter the vertical top concentration (uM):'),
        ('y_dilution', 'Enter the vertical dilution factor (x-fold):'),
    ]

    # Loop that creates labels and entry boxes for each field 
    # Fields: (key = internal identifier, text)
    for key, text in fields:
        ttk.Label(left_frame, text=text).pack() # left_frame is the main window, adds labels
        entry = ttk.Entry(left_frame) # adds entry widgets to left_frame
        entry.pack() # adds entries in vertical order
        entries[key] = entry # stores in dictionary

    # Image Selection for Genotype 1
    ttk.Label(left_frame, text="Select Genotype 1 Image Directory").pack()
    ttk.Button(left_frame, text="Select Genotype 1 Images", command=select_image_dir).pack()

    # Genotype 2
    ttk.Label(right_frame, text="Genotype 2 Parameters", font=('Arial', 12, 'bold')).pack()

    # Genotype 2 Fields
    entries_2 = {}
    fields_2 = [
        ('title_name_2', 'Enter the title of the experiment:'),
        ('x_drug_2', 'Enter the drug name for the horizontal gradient:'),
        ('x_top_conc_2', 'Enter the horizontal top concentration (uM):'),
        ('x_dilution_2', 'Enter the horizontal dilution factor (x-fold):'),
        ('y_drug_2', 'Enter the drug name for the vertical gradient:'),
        ('y_top_conc_2', 'Enter the vertical top concentration (uM):'),
        ('y_dilution_2', 'Enter the vertical dilution factor (x-fold):'),
    ]

    for key, text in fields_2:
        ttk.Label(right_frame, text=text).pack()
        entry = ttk.Entry(right_frame)
        entry.pack()
        entries_2[key] = entry

    # Image Selection for Genotype 2 
    ttk.Label(right_frame, text="Select Genotype 2 Image Directory (optional)").pack(pady=(10, 0))
    ttk.Button(right_frame, text="Select Genotype 2 Images", command=select_image_dir_2).pack()

    # This dictionary will hold the result
    gui_inputs = {}
    
    # Callback function to use when the form is submitted
    def submit():
        for key, entry in entries.items():
            gui_inputs[key] = entry.get()
        for key, entry in entries_2.items():
            gui_inputs[key] = entry.get()
        gui_inputs['image_dir'] = image_dir
        gui_inputs['image_dir_2'] = image_dir_2
        root.destroy()
    
    # Create button to submit form
    submit_button = ttk.Button(root, text='Submit', command=submit)
    submit_button.pack(pady = 20)
    
    # Runs GUI
    root.mainloop()
    logger.debug('GUI submitted.')
    return gui_inputs

# Assign the synergy_gui output to script variable
gui_inputs = synergy_gui()

# Assign GUI inputs to script variables
title_name = gui_inputs['title_name']
x_drug = gui_inputs['x_drug']
x_top_conc = float(gui_inputs['x_top_conc'])
x_dilution = float(gui_inputs['x_dilution'])
y_drug = gui_inputs['y_drug']
y_top_conc = float(gui_inputs['y_top_conc'])
y_dilution = float(gui_inputs['y_dilution'])

title_name_2 = gui_inputs['title_name_2']
x_drug_2 = gui_inputs['x_drug_2']
x_top_conc_2 = float(gui_inputs['x_top_conc_2'])
x_dilution_2 = float(gui_inputs['x_dilution_2'])
y_drug_2 = gui_inputs['y_drug_2']
y_top_conc_2 = float(gui_inputs['y_top_conc_2'])
y_dilution_2 = float(gui_inputs['y_dilution_2'])

# Two Image Directories for Dual Genotype Mode
image_dir = gui_inputs.get('image_dir')
image_dir_2 = gui_inputs.get('image_dir_2')

# Determine Mode: Dual or Single
if image_dir_2 and image_dir_2.strip() != "":
    mode = "dual"
else:
    mode = "single"

logger.info(f"Analysis mode detected: {mode}")
logger.debug('GUI inputs assigned to variables.')

# Calculate x and y concentration gradients
x_doses_1 = tb.dose_range_x(x_top_conc, x_dilution)
logger.debug('x_doses gradient calculated.')
y_doses_1 = tb.dose_range_y(y_top_conc, y_dilution)
logger.debug('y_doses gradient calculated.')

# Calculate x and y concentration gradients for genotype 2
x_doses_2 = tb.dose_range_x(x_top_conc_2, x_dilution_2)
logger.debug('x_doses_2 gradient calculated.')
y_doses_2 = tb.dose_range_y(y_top_conc_2, y_dilution_2)
logger.debug('y_doses_2 gradient calculated.')

def process_genotype(df_cp, x_doses, y_doses, title_name):
    # Load the CellProfiler Counts Into a DataFrame and Rename Wells
    df = df_cp.copy()
    df.drop('ImageNumber', axis=1, inplace=True)

    # Map TIFF File Names to Well IDs
    df.columns = ['nuclei', 'well']
    df['well'] = df['well'].apply(lambda x: tb.rename_wells(x, tb.wells))
    logger.debug('CellProfiler output rows renamed to well names.')

    # Extract Row/Column for Each Well for Vectorized Grouping
    df[['Row','Column']] = df['well'].str.extract(r'^([B-G])(\d+)$')
    logger.debug('Row and column labels extracted from well names.')

    # Group Triplicates and Calculate Statistics
    stats = (
        df
        .groupby('well')['nuclei']
        .agg(mean='mean', std='std')
        .reindex(tb.wells)
    )
    well_means = stats['mean'].tolist()
    well_std = stats['std'].tolist()
    logger.debug('Nuclei mean and standard deviation calculated.')

    # Map Concentrations/Doses
    tb_wells = [str(i) for i in range(2,12)]
    column_concentrations = dict(zip(
        tb_wells,
        [0] + x_doses
    ))
    row_concentrations = dict(zip(
        ['B','C','D','E','F','G'],
        [0] + y_doses
    ))

    # Map Concentrations to Rows and Columns
    rows = stats.index.str[0].tolist()
    columns = stats.index.str[1:].tolist()
    row_conc = [row_concentrations[r] for r in rows]
    col_conc = [column_concentrations[c] for c in columns]
    logger.debug('Row and column concentrations mapped.')

    # Normalize Well Means to Vehicle (B2) Mean
    vehicle = stats['mean'].iloc[0]
    normalized = stats['mean'] / vehicle

    # Frame the Data
    df_stats = pd.DataFrame({
        'Well': tb.wells,
        'Mean': stats['mean'],
        'Standard Deviation': stats['std'],
        'Normalized Mean': normalized,
        'Row Drug Concentration': row_conc,
        'Column Drug Concentration': col_conc
    })

    return df_stats, row_concentrations, column_concentrations

# Define File Path for synergy_output Subfolder
synergy_output_dir = base_dir / 'synergy_output'
synergy_output_dir.mkdir(exist_ok=True)
logger.debug('CellPyAbility/synergy_output/ identified or created and identified.')

# Run Genotype 1
df_cp_1, cp_csv_1 = tb.run_cellprofiler(image_dir)
df_stats_1, row_conc_1, col_conc_1 = process_genotype(
    df_cp_1, x_doses_1, y_doses_1, title_name)

# Save the Experiment Viability Statistics as a .csv
df_stats_1.to_csv(synergy_output_dir / f'{title_name}_synergy_stats.csv', index=False)
logger.info(f'{title_name} synergy stats saved to synergy_output')
tb.rename_counts(cp_csv_1, synergy_output_dir / f'{title_name}_synergy_counts.csv')

# Genotype 2
if mode == "dual":
    df_cp_2, cp_csv_2 = tb.run_cellprofiler(image_dir_2)

    df_stats_2, row_conc_2, col_conc_2 = process_genotype(
        df_cp_2, x_doses_2, y_doses_2, title_name_2
    )

    # Save the Experiment Viability Statistics as a .csv
    df_stats_2.to_csv(synergy_output_dir / f'{title_name_2}_synergy_stats.csv', index=False)
    logger.info(f'{title_name_2} synergy stats saved to synergy_output')
    tb.rename_counts(cp_csv_2, synergy_output_dir / f'{title_name_2}_synergy_counts.csv')

#WIP HERE

# # Initialize a list to store Bliss independence results
# bliss_results = []

# # Pull viability values from df_stats to calculate Bliss Independence
# for _, row in df_stats.iterrows():
#     well_name = row['Well']
#     observed_combined_effect = row['Normalized Mean']
#     logger.debug('Normalized means pulled from df_stats.')
    
#     # Determine the x-alone and y-alone effects based on the well name
#     if well_name[0] in 'BCDEFG' and well_name[1:] in '234567891011':
#         x_effect = df_stats.loc[df_stats['Well']=='B'+well_name[1:], 'Normalized Mean'].iloc[0]
#         y_effect = df_stats.loc[df_stats['Well']==well_name[0]+'2', 'Normalized Mean'].iloc[0]
#         logger.debug('x_effect and y_effect identified.')

#         # Calculate the expected combined effect
#         expected_combined_effect = x_effect * y_effect
#         logger.debug('Expected combined effects calculated.')

#         # Calculate the Bliss independence
#         bliss_independence = expected_combined_effect - observed_combined_effect
#         logger.debug('Bliss independence scores calculated.')

#         bliss_results.append({
#             'Well': well_name,
#             'Expected Combined Effect': expected_combined_effect,
#             'Observed Combined Effect': observed_combined_effect,
#             'Bliss Independence': bliss_independence
#         })

# # Convert the results to a DataFrame
# df_bliss = pd.DataFrame(bliss_results)

# # Add 'Row Drug Concentration' and 'Column Drug Concentration' to df_bliss
# df_bliss['Row Drug Concentration'] = df_bliss['Well'].str[0].map(row_concentrations)
# df_bliss['Column Drug Concentration'] = df_bliss['Well'].str[1:].map(column_concentrations)

# # Create pivot tables for plotting and saving
# normalized_means_pivot = df_stats.pivot(
#     index='Row Drug Concentration', columns='Column Drug Concentration', values='Normalized Mean'
# )
# bliss_independence_pivot = df_bliss.pivot(
#     index='Row Drug Concentration', columns='Column Drug Concentration', values='Bliss Independence'
# )

# # Convert pivot tables to numpy arrays
# cell_survival = normalized_means_pivot.values
# bliss_independence = bliss_independence_pivot.values

# # Save viability and Bliss matrices as .csv files
# normalized_means_pivot.to_csv(synergy_output_dir / f'{title_name}_synergy_ViabilityMatrix.csv')
# logger.info(f'{title_name} viability matrix saved to synergy_output.')
# bliss_independence_pivot.to_csv(synergy_output_dir / f'{title_name}_synergy_BlissMatrix.csv')
# logger.info(f'{title_name} Bliss score matrix saved to synergy_output.')

# # Extract x and y values from the pivot tables
# x_values = normalized_means_pivot.columns.values
# y_values = normalized_means_pivot.index.values

# # Replace zero values in x_values and y_values with a small positive number of fixed logarithmic distance
# min_x_value = x_values[x_values>0].min()
# min_y_value = y_values[y_values>0].min()

# x_values = np.where(x_values==0, min_x_value/(x_values[3]/x_values[2]), x_values)
# y_values = np.where(y_values==0, min_y_value/(y_values[3]/y_values[2]), y_values)

# x_tickvals = np.unique(np.concatenate(([min_x_value/(x_values[3]/x_values[2])], x_values)))
# y_tickvals = np.unique(np.concatenate(([min_y_value/(y_values[3]/y_values[2])], y_values)))

# x_ticktext = ['0'] + [f'{val:.1e}' for val in x_tickvals[1:]]
# y_ticktext = ['0'] + [f'{val:.1e}' for val in y_tickvals[1:]]
# logger.debug('Replaced zero with non-zero values a fixed log distance from the minimum concentrations.')

# # Create the 3D surface plot
# fig = go.Figure(data=[go.Surface(z=cell_survival, x=x_values, y=y_values, surfacecolor=bliss_independence, colorscale='jet_r', cmin=-0.3, cmax=0.3, colorbar=dict(title='Bliss Independence'))])
# logger.debug('3D surface plot created.')

# # Update layout to set x and y axes to logarithmic scale
# fig.update_layout(title=str(title_name), scene=dict(xaxis=dict(title=x_drug, type='log', ticktext=x_ticktext, tickvals=x_tickvals), yaxis=dict(title=y_drug, type='log', ticktext=y_ticktext, tickvals=y_tickvals), zaxis=dict(title='Relative Cell Survival', range=[0,np.max(cell_survival)])))
# logger.debug('x and y axes changed to log, experiment details added as labels.')

# # Rename CellProfiler output with experiment title and save in synergy_output
# counts_csv = synergy_output_dir / f'{title_name}_synergy_counts.csv'

# tb.rename_counts(cp_csv, counts_csv)
# logger.info(f'{title_name} raw counts saved to synergy_output.')

# # Save the interactable plot as an HTML
# fig.write_html(synergy_output_dir / f'{title_name}_synergy_plot.html')
# logger.info(f'{title_name} synergy plot saved to synergy_output.')

# fig.show()