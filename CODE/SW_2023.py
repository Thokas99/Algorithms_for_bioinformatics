# title: "Algorithms for Bioinformatics Project 2023, QCB, Unitn"
# version : 2.4v Multi-result
# author: "Thomas Sirchi"
# date: "26/05/2023"

#Packages essential to the code:
#pip install questionary
import questionary
#pip install click
import click
#pip install textwrap
import textwrap
from collections import defaultdict
#pip install itertools
import itertools
#pip install numpy
import numpy as np
#pip install matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
#pip install seaborn
import seaborn as sns
#pip install tabulate
from tabulate import tabulate
#pip install pandas
import pandas as pd
#pip install biopython
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight


##################################################################################################################
##################################################################################################################

#### IMPUT CHECK ####


def is_valid_dna_sequence(sequence):
    """
    Check if a given DNA sequence contains only valid nucleotides (A, T, G, C).

    Args:
     - sequence (str): the DNA sequence to be validated.

    Returns:
     - bool: True if the sequence contains only valid nucleotides, False otherwise.
    """

    # Define a set of valid nucleotide characters.
    valid_chars = set('ATGC')

    # Convert the input sequence to a string  and uppercase, in case it's not already a string.
    sequence = str(sequence).upper()

    # Loop through each character in the sequence.
    for char in sequence:
        # Check if the character is not a valid nucleotide.
        if char not in valid_chars:
            # Print a warning message indicating which character is invalid.
            print(f"\n{char}\t <-- is not a nucleotide \n")

            # Set a flag to indicate that the sequence is not valid.
            is_good = False
            # Exit the loop.
            break
    else:
        # If the loop completes without breaking, then all characters are valid nucleotides.
        # Set the flag to indicate that the sequence is valid.
        is_good = True

    # Return the flag indicating whether the sequence is valid or not.
    return is_good


def Input_check(seq_a,seq_b):
    """
    This function checks if the input sequences are valid DNA sequences, and prompts the user to enter new sequences if they are invalid or empty.

    Args:
     - seq_a (str): the first DNA sequence to be checked
     - seq_b (str): the second DNA sequence to be checked

    Returns:
     - seq_a (str): the validated first DNA sequence
     - seq_b (str): the validated second DNA sequence
    """
    while True:
        # Check if either sequence is empty
        if seq_a == "":
            print("Sequence A must not be empty.\n")
            # Prompt user to enter new sequences or choose from a list of demos
            seq_a = questionary.autocomplete("Enter a new sequence A or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
        # Check if either sequence is empty
        elif seq_b == "":
            print("Sequence B must not be empty.\n")
            seq_b = questionary.autocomplete("Enter a new sequence B or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
        
        # Check if either sequence is not a valid DNA sequence
        elif not is_valid_dna_sequence(seq_a) and not is_valid_dna_sequence(seq_b):
            print("Neither the first and the second sequences are a valid DNA sequences.\n")
            # Prompt user to enter new sequences or choose from a list of demos
            seq_a = questionary.autocomplete("Enter a new sequence A or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
            seq_b = questionary.autocomplete("Enter a new sequence B or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
        
        # Check if the first sequence is not a valid DNA sequence
        elif not is_valid_dna_sequence(seq_a) and is_valid_dna_sequence(seq_b):
            print("The first sequence is not a valid DNA sequence.\n")
            # Prompt user to enter a new sequence a or choose from a list of demos
            seq_a = questionary.autocomplete("Enter a new sequence A or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
        
        # Check if the second sequence is not a valid DNA sequence
        elif not is_valid_dna_sequence(seq_b) and is_valid_dna_sequence(seq_a):
            print("The second sequence is not a valid DNA sequence.\n")
            # Prompt user to enter a new sequence b or choose from a list of demos
            seq_b = questionary.autocomplete("Enter a new sequence B or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()

        # If both sequences are valid, break out of the loop and return the sequences
        else:
            break
    return seq_a, seq_b


def ask_to_repeat(seq_a, seq_b):
    """
    Asks the user whether to run the Smith-Waterman function again, and whether to use new sequences.

    args:
     - seq_a (str): The first sequence to be aligned.
     - seq_b (str): The second sequence to be aligned.
    """

    # Continuously loop until the user chooses to exit
    while True:
        # Ask the user whether to run the function again
        answer = questionary.select("Do you want to run the function again?", choices=["Yes", "No"]).ask()

        # If the user chooses to exit, break out of the loop
        if answer == 'No':
            print("Thank you for using the Smith-Waterman algorithm script!!!","\U0001F601")
            break
        # If the user chooses to run the function again, ask whether to use new sequences
        elif answer == 'Yes':
            answer2 = questionary.select("Do you want to use new sequences or parameters?", choices=["Yes, just new sequences",
                                                                                                     "Yes, just new parameters",
                                                                                                     "Yes, new sequences and parameters", 
                                                                                                     "No"]).ask()

            # If the user chooses to use the same sequences, run the function again with the same inputs
            if answer2 == 'No':
                smith_waterman(seq_a, seq_b,user_input = False)
                break

            # If the user chooses to use new sequences or parameters, ask for the corresponding and run the function again with the new inputs
            elif answer2 == 'Yes, just new sequences':
                new_seq_a = questionary.autocomplete("Enter a new sequence A or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
                new_seq_b = questionary.autocomplete("Enter a new sequence B or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
                smith_waterman(new_seq_a, new_seq_b,user_input = False)
                break

            elif answer2 == 'Yes, just new parameters':
                smith_waterman(seq_a, seq_b,user_input = False,params=True)
                break

            elif answer2 == 'Yes, new sequences and parameters':
                new_seq_a = questionary.autocomplete("Enter a new sequence A or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
                new_seq_b = questionary.autocomplete("Enter a new sequence B or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
                smith_waterman(new_seq_a, new_seq_b,user_input = False,params=True)
                break

##################################################################################################################
##################################################################################################################

#### PRINT ####

def input_stats(seq_a, seq_b):
    """
    Calculates various statistics and displays information about the input sequences.

    Args:
        seq_a (str): The first DNA sequence.
        seq_b (str): The second DNA sequence.
    """

    # Convert sequences to uppercase and create Seq objects
    seq_a = Seq(seq_a.upper())
    seq_b = Seq(seq_b.upper())

    # Calculate GC fraction for each sequence (0-1)
    gc_a = round(gc_fraction(seq_a), 2)
    gc_b = round(gc_fraction(seq_b), 2)

    # Calculate weight of the sequences
    weight_a = molecular_weight(seq_a, "DNA")
    weight_b = molecular_weight(seq_b, "DNA")

    # Calculate Hamming distance between sequences
    hamming_distance = sum(1 for a, b in zip(seq_a, seq_b) if a != b)

    # Create a dictionary of nucleotide counts for each sequence
    count_a = dict(zip(["A", "T", "C", "G"], [seq_a.count("A"), seq_a.count("T"), seq_a.count("C"), seq_a.count("G")]))
    count_b = dict(zip(["A", "T", "C", "G"], [seq_b.count("A"), seq_b.count("T"), seq_b.count("C"), seq_b.count("G")]))

    # Prepare statistics for display
    stats = [
        ['GC content of A', count_a["G"] + count_a["C"], f"{gc_a}"],
        ['GC content of B', count_b["G"] + count_b["C"], f"{gc_b}"]]

    stats_2 =[
        ["Weight of A", f"{weight_a:.2f}g/mol"],
        ["Weight of B", f"{weight_b:.2f}g/mol"],
        ['Hamming distance', hamming_distance,]]
   

    # Print statistics using tabulate for a formatted table
    print("\nStatistics on the input sequences:")
    print(tabulate(stats, numalign="center", tablefmt="simple_grid", headers=['Statistic', 'Count', 'GC fraction (0-1)']), "\n")
    print(tabulate(stats_2, numalign="center", tablefmt="simple_grid", headers=['Statistic', 'Count']), "\n")

    # Create a pandas DataFrame to display nucleotide counts for each sequence
    data = {"Seq A": count_a, "Seq B": count_b}
    df = pd.DataFrame(data).transpose()

    # Display nucleotide counts using tabulate for a formatted table
    print("Nucleotide counts:\n", tabulate(df, numalign="center", headers=["Sequence","A", "T", "C", "G"], tablefmt="grid"), '\n')


def print_matrix(matrix ,positions, seq_a, seq_b, threshold = 30, heatmap = True):
    """Prints a matrix up to a certain threshold value.
    If the matrix is larger than the threshold, it only prints the first rows
    and provides an aesthetic message to indicate that the matrix is too long.

    Args:
     - seq_a (str): Aligned sequence A.
     - seq_b (str): Aligned sequence B.
     - matrix (numpy.ndarray): The matrix to be printed.
     - threshold (int): The maximum number of rows to be printed.
     - heatmap (bool): If True, plot a heatmap of the matrix using seaborn.
     - positions (list): The list of posiion visited in the matrix to create the arrows

    """
    # Getting the sequence ready to be plotted on the heatmap
    plot_a = " "+ seq_a
    plot_b = " "+ seq_b
    num_rows = matrix.shape[0]
    if num_rows <= threshold:
        # Print the entire matrix if it's smaller than the threshold
        print("\nMatrix:\n" )
        # Create a DataFrame from the matrix data and axis labels to show the sequance
        print(tabulate(pd.DataFrame(matrix, index=list(plot_a)), tablefmt="simple_grid",numalign="center", headers = list(plot_b)),"\n")
        if heatmap:
            # Create list of 148 colors to use and a count to fifferenciate each arrow
            colors = list(mcolors.CSS4_COLORS.values())
            # List to keep track of the arrows in the legend
            arrow_legends = []
            # Plot a heatmap
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.heatmap(matrix, annot=True, cmap='coolwarm', annot_kws={'fontsize':12}, 
                        square=True, linewidths=0.01, cbar_kws={"shrink": 0.8}, cbar= False,
                        xticklabels=["{}".format(n) for n in plot_b],
                        yticklabels=["{}".format(n) for n in plot_a])
            ax.xaxis.tick_top()
            plt.xlabel('Sequence B')
            plt.ylabel('Sequence A')
            plt.title('Smith-Waterman Scoring Matrix')
            fig.savefig('my_plot.png', dpi=300, bbox_inches='tight')
            # Add the rows to show each path in the heathmap 
            # if condition to avoid index out of posistion for more than 148 arrows
            if len(positions) <= 100:
                for i in range(0,len(positions)+1):
                    n = np.random.randint(100)
                    my_col = colors[n]
                    for pos in range(len(positions[i])-1):
                        start = positions[i][pos]
                        end = positions[i][pos+1]
                        ax.annotate('', xy=end, xytext=start, arrowprops=dict(arrowstyle='->', lw=2, color=my_col))
                    if i > 0: # Starting from 1 because i is used to iterate a dictionary and not a list
                        # Create a Line2D object for the arrow legend
                        arrow_legend = plt.Line2D([0], [0], marker='o', color=my_col, label=str(i) + "° result", markerfacecolor= my_col, markersize=5)
                        # Append the arrow legend to the list
                        arrow_legends.append(arrow_legend)
            else:
                for i in range(0,len(positions)+1):
                    for pos in range(len(positions[i])-1):
                        start = positions[i][pos]
                        end = positions[i][pos+1]
                        ax.annotate('', xy=end, xytext=start, arrowprops=dict(arrowstyle='->', lw=2, color="r"))
                    
                # Create a Line2D object for the arrow legend
                arrow_legend = plt.Line2D([0], [0], marker='o', color="r", label="More than 100 results!", markerfacecolor= "r", markersize=5)
                # Append the arrow legend to the list
                arrow_legends.append(arrow_legend)

            # Add the arrow legends to the plot
            ax.legend(title='Alignments:', loc='center left', handles=arrow_legends, bbox_to_anchor=(1, 0.5))
            plt.show()
            fig.savefig('my_plot_arrows.png', dpi=300, bbox_inches='tight')
            plt.close()
    else:
        # Print only the first threshold rows
        print("\n",tabulate(matrix[:threshold][:threshold]),"\n")
        # Print an aesthetic message to indicate that the matrix is too longg
        print("... matrix is too long, only showing the first {} rows and columns".format(threshold))


def generate_alignment_string(aligned_seq_a, aligned_seq_b): 
    """
    Generates a string representing the alignment of two sequences.

    Args:
     - aligned_seq_a (str): The first aligned sequence.
     - aligned_seq_b (str): The second aligned sequence.

    Returns:
     - alignment_string (str): The string representing the alignment.
     - idents (int): The number of matching characters in the alignment.
     - gaps (int): The number of gap characters in the alignment.
     - mismatches (int): The number of mismatched characters in the alignment.
    """
    # Initialize variables to track the number of matches, gaps, and mismatches
    idents, gaps, mismatches = 0, 0, 0
    # Initialize an empty list to build the alignment string
    alignment_string = []
    
    # Loop through both sequences simultaneously
    for base_a, base_b in zip(aligned_seq_a, aligned_seq_b): 
        if base_a == base_b:
            # If the two characters match, add a '|' character to the alignment string
            alignment_string.append('|')
            idents += 1
        elif '-' in (base_a, base_b):
            # If either character is a gap, add a ' ' character to the alignment string
            alignment_string.append(' ')
            gaps += 1
        else:
            # If the characters do not match and are not gaps, add an 'X' character to the alignment string
            alignment_string.append('X')
            mismatches += 1
    # Combine the characters in the alignment string list into a single string
    visual_alig = ''.join(alignment_string)
    
    return visual_alig, idents, gaps, mismatches


def pretty_alignment(aligned_seq_a, aligned_seq_b, alignment_score,thresholds):
    """
    Print a pretty alignment of two sequences and related statistics.
    
    Args:
     - aligned_seq_a (str): Aligned sequence A.
     - aligned_seq_b (str): Aligned sequence B.
     - alignment_score (list): List of alignment scores
     - thresholds (dict) : Dict of thresholds to filter sequences 
    """

    # Get alignment statistics
    alignment_str, idents, gaps, mismatches = generate_alignment_string(aligned_seq_a, aligned_seq_b)



    if gaps >= thresholds["gaps"]:
        consecutive = "|"*thresholds["matches"]
        if  consecutive in alignment_str:
            # If no alignment is available, print a message and return
            if len(aligned_seq_a) == 0: 
                return print("There is no alignment available between the sequences")
            elif len(aligned_seq_a) >= 0 and len(aligned_seq_b) >= 0:
                ####  Print some stats ####
                alength = len(aligned_seq_a)
                # Create a list of lists containing the data
                stats = [
                        ['Identities', f"{idents}/{alength}", f"{idents / alength:.1%}"],
                        ['Gaps', f"{gaps}/{alength}", f"{gaps / alength:.1%}"],
                        ['Mismatches', f"{mismatches}/{alength}", f"{mismatches / alength:.1%}"],
                        ['Total score', alignment_score, '']
                        ]
                # Use the tabulate function to format the data and print it
                print(tabulate(stats,tablefmt= "simple_grid",headers=['Statistic', 'Count', 'Percentage']),"\n")

                # Determine the maximum sequence length
                max_length = max(len(aligned_seq_a), len(alignment_str), len(aligned_seq_b))

                # Define the chunk width as the maximum sequence length plus 10 characters for padding
                chunk_width = max_length + 10

                # Wrap aligned sequences into chunks
                seqA_chunks = textwrap.wrap(aligned_seq_a, width=chunk_width)
                align_chunks = textwrap.wrap(alignment_str, width=chunk_width)
                seqB_chunks = textwrap.wrap(aligned_seq_b, width=chunk_width)

                # Print sequences in chunks
                for i, (seqA_chunk, align_chunk, seqB_chunk) in enumerate(zip(seqA_chunks, align_chunks, seqB_chunks), start=1):
                    print(f"Seq_A {i*chunk_width-chunk_width+1:<4} {seqA_chunk}")
                    print(f"           {align_chunk}")
                    print(f"Seq_B {i*chunk_width-chunk_width+1:<4} {seqB_chunk}")
                    print()
        else:
            print("Result not conforming with the thresholds")
    else:
        print("Result not conforming with the thresholds")
##################################################################################################################
##################################################################################################################

#### SMITH_WATERMAN ####


def matrix(seq_a, seq_b, match_score=3, gap_cost=2,mismatch = -3):
    """
    Create an alignment matrix for two input sequences.

    Args:
     - seq_a (str): the first sequence to be aligned.
     - seq_b (str): the second sequence to be aligned.
     - match_score (int): the score for a match between two nucleotides (default is 3).
     - gap_cost (int): the penalty score for inserting a gap in the alignment (default is 2).
     - mismatch (int): the penalty score for a mismatch in the aligment (default is -3)

    Returns:
     - np.ndarray: a 2D numpy array containing the alignment matrix.
    """

    # Create an alignment matrix filled with zeros, with dimensions (len(seq_a) + 1) x (len(seq_b) + 1)
    alignment_matrix = np.zeros((len(seq_a) + 1, len(seq_b) + 1), int)

    # Loop over all indices of the alignment matrix (except the top row and left column)
    for i, j in itertools.product(range(1, alignment_matrix.shape[0]), range(1, alignment_matrix.shape[1])):
        # Calculate the scores for three possible moves (match/mismatch, gap in seq_a, gap in seq_b)
        match = alignment_matrix[i - 1, j - 1] + (match_score if seq_a[i - 1] == seq_b[j - 1] else mismatch)
        delete = alignment_matrix[i - 1, j] - gap_cost
        insert = alignment_matrix[i, j - 1] - gap_cost
        # Set the score at the current index of the alignment matrix to the maximum of the three possible moves or 0
        alignment_matrix[i, j] = max(match, delete, insert, 0)

    # Return the completed alignment matrix
    return alignment_matrix


def traceback(alignment_matrix, seq_a, seq_b,thresholds, match_score=3, gap_cost=2, mismatch =-3):
    """
    Given an alignment matrix and two sequences, performs traceback to find all possible alignments 
    with the highest score(s). Returns the aligned sequences, score, movement, and positions.

    Args:
     - alignment_matrix (numpy.ndarray): The alignment matrix.
     - seq_a (str): The first sequence.
     - seq_b (str): The second sequence.
     - thresholds (dict) : Dict of thresholds to filter sequences 
     - match_score (int): the score for a match between two nucleotides (default is 3).
     - gap_cost (int): the penalty score for inserting a gap in the alignment (default is 2).
     - mismatch (int): the penalty score for a mismatch in the aligment (default is -3)

    Returns:
     - all_aligned_seq_a (defaultdict): A dictionary of all aligned sequences for seq_a.
     - all_aligned_seq_b (defaultdict): A dictionary of all aligned sequences for seq_b.
     - plot_positions (defaultdict): A dictionary of positions for plotting.
     - max_scores (defaultdict): A dictionary of maximum scores for each alignment.
     - trace_len (defaultdict): A dictionary of len of the sequances for each alignment.
    """
    
    # get the index of the maximum value in the matrix
    max_value = abs((np.amax(alignment_matrix)/100)*thresholds["aligment_score_low"])
    print(f"\nStarting from the value {max_value}\n")
    max_value_index = np.argwhere(alignment_matrix >= max_value)
    print("The values in the matrix are at:\n", tabulate(max_value_index,numalign="center",tablefmt= "grid",headers= (" X "," Y ")),"\n")
    
    # initialize number of itarations (tmp)
    tmp = 0
    #inizialize the dictionary to store multiple results
    all_aligned_seq_a = defaultdict()
    all_aligned_seq_b = defaultdict()
    plot_positions = defaultdict(list)
    max_scores = defaultdict()
    trace_len = defaultdict()
    # trace the alignment path
    for tup in max_value_index:
        p_positions = list()
        # initialize empty aligned sequences
        aligned_seq_a, aligned_seq_b = '', ''
        tmp +=1
        i, j = tup[0], tup[1]
        score = alignment_matrix[i, j]
        while alignment_matrix[i, j] != 0:
            # check if the current score is from a match/mismatch
            if alignment_matrix[i, j] == alignment_matrix[i-1, j-1] + (match_score if seq_a[i-1] == seq_b[j-1] else mismatch):
                # add aligned characters and move diagonally up-left
                p_positions.append((j+0.5,i+0.5))
                aligned_seq_a = seq_a[i-1] + aligned_seq_a
                aligned_seq_b = seq_b[j-1] + aligned_seq_b
                i, j = i-1, j-1
            # check if the current score is from a gap in seq_a
            elif alignment_matrix[i, j] == alignment_matrix[i-1, j] - gap_cost:
                # add gap and move up
                aligned_seq_a = seq_a[i-1] + aligned_seq_a
                aligned_seq_b = '-' + aligned_seq_b
                i -= 1
                p_positions.append((j+0.5,i+1.5))
            # check if the current score is from a gap in seq_b
            else:
                # add gap and move left
                aligned_seq_a = '-' + aligned_seq_a
                aligned_seq_b = seq_b[j-1] + aligned_seq_b
                j -= 1
                p_positions.append((j+1.5,i+0.5))

            # add the current score to the score list

        # Add al the traceback to the respective dict 
        plot_positions[tmp] = p_positions
        all_aligned_seq_a[tmp] = aligned_seq_a
        all_aligned_seq_b[tmp] = aligned_seq_b
        max_scores[tmp] = score
        trace_len[tmp] = len(aligned_seq_a) #both aligned_seq_a and aligned_seq_b, so I can only register aligned_seq_a

    # Convert dictionaries to dataframes
    df_seq_a = pd.DataFrame.from_dict(all_aligned_seq_a, orient='index', columns=['Aligned Seq A'])
    df_seq_b = pd.DataFrame.from_dict(all_aligned_seq_b, orient='index', columns=['Aligned Seq B'])

    # Merge dataframes on their indices
    df_merged = pd.merge(df_seq_a, df_seq_b, left_index=True, right_index=True)

    # Print the merged dataframe
    print("All the alignments:\n")
    print(df_merged)
    # return the aligned sequences, score, movement, and positions
    return all_aligned_seq_a, all_aligned_seq_b, plot_positions,max_scores,trace_len


def smith_waterman(seq_a="", seq_b="", match_score=3, gap_cost=2, mismatch = -3,user_input = True, ask = True,params = False, graphic=True):
    """
    The function smith_waterman performs the Smith-Waterman algorithm for local sequence alignment between two input sequences. 
    The algorithm generates an alignment matrix based on the match_score and gap_cost parameters, and then performs a traceback to find the optimal alignment.

    Args:
     - seq_a: the first input sequence
     - seq_b: the second input sequence
     - match_score (int): the score for a match between two nucleotides (default is 3).
     - gap_cost (int): the penalty score for inserting a gap in the alignment (default is 2).
     - mismatch (int): the penalty score for a mismatch in the aligment (default is -3)
     - user_input: a boolean flag indicating whether to prompt the user to ask the sequences (default=True)
     - params: a boolean flag indicating whether to prompt the user to ask if he wants to modify the parameters (default=True)
     - ask: a boolean flag indicating whether to prompt the user to repeat the function (default=True)
     - graphic: a boolean flag indicating whether to display a graphical representation of the alignment matrix (default=True)
     """
    
    # Set seed
    np.random.seed(710)
    # Set thresholds for the modification required by the exam
    thresholds = {"gaps":2,"matches": 4,"aligment_score_low":50}

    #ask for the sequances
    if user_input:
        print(tabulate([["Welcome to the Smith-Waterman algorithm!"], 
              ["This algorithm is named after Smith and Waterman, who developed it."],
              ["Thomas Sirchi has written the code for this script."]], tablefmt="simple_grid"))
        seq_a = questionary.autocomplete("Enter a sequence A or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
        seq_b = questionary.autocomplete("Enter a sequence B or chose a demo:", choices=["ATCGCTAGCG","AGTTCGCTGA","GCTAGCTAGC","ATGCTAGCTA","CTAGCTGAGC","TGTTACGGAG","GGTTGACTAT"]).ask()
        # Ask the user if he whats to change the param
        answer = questionary.select("Do you want to input the parameters and thresholds for the algoritm? \nDefaults parameters:[ match_score=3, gap_cost=2, mismatch = -3]\nDefaults thresholds:[ total gaps = 2, consegutive matches = 4, Lowest aligment score = 50/100 of maximum aligment score ]", choices=["Yes", "No"]).ask()

        # If the user chooses to exit, break out of the loop
        if answer == 'No':
            print("Default parameters will be applied")
            user_input = False
        # If the user chooses to run the function again, ask whether to use new sequences
        elif answer == 'Yes':
            params = True
            answer = questionary.select("Do you want to apply the thresholds for the algoritm? Defaults:[ gaps = 2, consegutive matches = 4 ]", choices=["Yes, custom", "No, defaults"]).ask()
            # If the user chooses to exit, break out
            if answer == 'No, defaults':
                print("Default thresholds will be applied")
            # if the user wants let the inputs 
            if answer == "Yes, custom":
                gaps = abs(click.prompt("Please insert a positive gap number", type=int)) # must be positive 
                matches = abs(click.prompt("Please insert a positive number of consegutive matches ", type=int)) # must be positive
                aligment_score_low = abs(click.prompt("Please insert a positive threshold percentage to calculate the lowest score", type=int)) # must be positive
                thresholds = {"gaps":gaps,"matches": matches,"aligment_score_low": aligment_score_low}
                user_input = False

    if params:   
        match_score = abs(click.prompt("Please insert a positive match score", type=int)) # must be positive 
        gap_cost = abs(click.prompt("Please insert a positivie gap cost", type=int)) # must be positive 
        mismatch = -abs(click.prompt("Please insert a negative mismatch cost", type=int)) # must be negative 
        params=False

    # Check the input
    true_seq_a, true_seq_b = Input_check(seq_a, seq_b)

    # Convert input sequences to uppercase
    true_seq_a, true_seq_b = true_seq_a.upper(), true_seq_b.upper()

    # Print stats
    input_stats(true_seq_a,true_seq_b)

    # Calculate the alignment matrix
    alignment_matrix = matrix(true_seq_a, true_seq_b, match_score, gap_cost, mismatch)

    # Traceback to find the optimal alignment
    aligned_seq_a, aligned_seq_b, p_positions,max_score,trace_len = traceback(alignment_matrix, true_seq_a, true_seq_b,thresholds,match_score, gap_cost, mismatch)

    # Sort the aligments by len as requested by the exam 
    sorted_trace_len = dict(sorted(trace_len.items(), key=lambda x: x[1], reverse=True))
    
    for i in list(sorted_trace_len.keys()):
        print(f"\n{i}° aligment possible between the sequences ####\n")
        pretty_alignment(aligned_seq_a[i],aligned_seq_b[i],max_score[i],thresholds)

    # Print the alignment matrix if graphic=True to help visualize the alignment
    if graphic:
        print_matrix(alignment_matrix, p_positions, true_seq_a, true_seq_b)

    # Ask the user whether to terminate or run the function again
    if ask:
        ask_to_repeat(true_seq_a, true_seq_b)
    
##################################################################################################################   
##################################################################################################################

smith_waterman()
