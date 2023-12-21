##################################################################################################################
##################################################################################################################

#### FASTA IMPORT (BETA) ####

def analyze_fasta_file(file_path, show_plots=True):
    # Packages for the FASTA import (If necessary)
    import pandas as pd
    from Bio import SeqIO,SeqUtils
    from Bio.Seq import Seq
    from tabulate import tabulate
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Open the FASTA file and parse the records
    records = list(SeqIO.parse(file_path, "fasta"))

    # Extract the sequence strings
    sequence = "".join([str(record.seq) for record in records])

    # Extract statistics from the records
    seq = Seq(sequence)
    num_records = len(records)
    lengths = [len(record.seq) for record in records]
    total_length = sum(lengths)
    avg_length = total_length / num_records
    min_length = min(lengths)
    max_length = max(lengths)

    # Generate and show some plots
    if show_plots:
        """# Plot a histogram of sequence lengths
        sns.histplot(lengths, bins=50)
        plt.xlabel("Sequence Length")
        plt.ylabel("Frequency")
        plt.title("Sequence Length Distribution")
        plt.show()"""

        # Plot a bar chart of nucleotide frequencies in the first sequence
        seq_string = str(records[0].seq)
        base_counts = dict.fromkeys(["A", "C", "G", "T"], 0)
        for base in seq_string:
            if base in base_counts:
                base_counts[base] += 1
        sns.barplot(x=list(base_counts.keys()), y=list(base_counts.values()))
        plt.xlabel("Nucleotide")
        plt.ylabel("Frequency")
        plt.title("Nucleotide Frequencies in First Sequence")
        plt.show()

    # Create a pandas DataFrame to store the statistics
    stats_df = pd.DataFrame({
        "ID": ["Total"] + [record.id for record in records],
        "Num Records": [num_records] + [1] * num_records,
        "Total Length": [total_length] + lengths,
        "A content": base_counts["A"],
        "C content": base_counts["C"],
        "G content": base_counts["G"],
        "T content": base_counts["T"],
        "GC content": "{:.2f}%".format(SeqUtils.GC(seq)),
        "Avg Length": [avg_length] + [len(record.seq) for record in records],
        "Min Length": [min_length] + [len(record.seq) for record in records],
        "Max Length": [max_length] + [len(record.seq) for record in records]
        
    })

    # Print the summary statistics
    print(tabulate(stats_df, headers='keys', tablefmt='psql'))


    # Return the sequence strings
    return sequence