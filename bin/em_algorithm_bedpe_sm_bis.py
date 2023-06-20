#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import os
import errno
import multiprocessing as mp
import time
import sys
def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception

# This code implements an Expectation-Maximization (EM) algorithm with a Bayesian prior for distributing multi-mapping reads (i.e., sequencing reads that map to multiple locations in the genome) to specific targets (e.g., genes). 
# The purpose of this algorithm is to assign a weight to each target using single mapping to targets as initial weight and optimizing it through prior/posterior probabilities of multi-mapping reads.
# The code is expected to read an input file in a specific format, preprocess the data, compute the prior probability of mapping to each target using single-mappers information, and then iteratively update the prior probabilities using the posterior probabilities computed from the multi-mappers data until convergence.
# It appears that the code is implementing the necessary functions to perform these tasks. 

# Expectation Maximization algorithm with Bayesian prior for multi-mapping reads distribution to targets:

# The input is a BEDPE file with both single and multi mappers

# Parse arguments
# Input file: a BEDPE file without header. Columns are:
# 1) read_full_name
# 2) read_id
# 3) target_id
# Max iterations of the EM algorithm (default = 1000)
# Convergence value of the EM (default 1e-6)
# Maximum iterations of the EM algorithm (default = 1000)
# Number of processors (default = 4)

# To run the script tipe python EM.py -i input_file -o output_file -m number_of_processes -c convergence_threshold i.e the convergenceilon value

def parse_args():
    parser = argparse.ArgumentParser(description='Expectation Maximization algorithm with Bayesian prior for multi-mapping reads distribution to targets')
    parser.add_argument('-i', '--input_file', help='Input file: a BEDPE file without header. Columns are: 1) chromosome of the first read, 2) start position of the first read, 3) end position of the first read, 4) chromosome of the second read, 5) start position of the second read, 6) end position of the second read, 7) read ID combined with the Hi:i: tag to distinguish individual multi-mapping reads. Original read ID is extracted by removing \'_HI:.*\' from the read ID, 8) Mapping quality of the first read, 8) strand of the first read, 9) strand of the second read, 10) target ID', required=True)
    parser.add_argument('-o', '--output_file', help='Output to the folder where the output file will be saved', required=True)
    parser.add_argument('-m', '--max_iter', help='Max iterations of the EM algorithm (default = 1000)', default=1000, type=int)
    parser.add_argument('-c', '--convergence', help='Convergence value of the EM (default 1e-6)', default=1e-6, type=float)
    args = parser.parse_args()
    return args

# read_input_file function to read the input file and preprocess the data
# Given a path to a BEDPE file, it returns one dataframe where we keep multimappers with taget_id covered by single mappers.
# After the inport the output df will have number of columns = 13 and number of rows = number of multimappers with taget_id covered by single mappers.
def read_input_file(path_bedpe):
    # path_bedpe is the path to the BEDPE file
    # Read the input file
    df = pd.read_csv(path_bedpe, sep='\t', header=None)
    # Check that the input file has 3 columns rread_full , ead_id and target_id 
    if df.shape[1] != 3:
        raise ValueError('The input file should have 3 columns')
    # Rename all columns with the vactor of identifiers
    df = df.rename(columns={0: 'read_full_name', 1: 'read_id', 2: 'target_id'}) 
    # Count the number of times each read_id appears in the dataframe
    df['counts_per_read_id'] = df.groupby('read_id')['read_id'].transform('count')
    # Get single mappers
    df_single = df[df['counts_per_read_id'] == 1]
    # Get multi mappers
    df_multi = df[df['counts_per_read_id'] > 1]
    # Filter multi mappers to keep only those with target_id covered by single mappers
    df_multi = df_multi[df_multi['target_id'].isin(df_single['target_id'])]
    # Filter single mappers to keep only those with target_id covered by multi mappers
    df_single = df_single[df_single['target_id'].isin(df_multi['target_id'])]
    # Join the two dataframes back together
    df = pd.concat([df_single, df_multi])
    # Sort the dataframe by read_id and Hi
    df = df.sort_values(by=['read_id', 'counts_per_read_id'])
    # Reset the index
    df = df.reset_index(drop=True)
    # Return the dataframe
    return df

# We aim at updating the prior probability associated to each target given the whole set of multi-mappers and the current prior probability of each target.
# Considering the Bayesian theory i.e. Posterior P(tj|r_all) = P(r_all|tj) * P(tj) / P(r_all), we have:
# 1. P(r_all|tj) = P(r1|tj) * P(r2|tj) * ... * P(rn|tj) and is named as the likelihood of the reads given the target tj.
# 2. P(tj) = prior probability of the target tj
# 3. P(r_all) = P(r1) * P(r2) * ... * P(rn) and is named as the likelihood of the reads which is a constant and can be left out.
# We basically resolve : P(tj|r_all) = P(r1|tj) * P(r2|tj) * ... * P(rn|tj) * P(tj) as a quantity that at each iteration will be normalized to 1 and used as the new prior probability of the target tj.

# # Consider the following example in which we have 4 targets (t1,t2,t3,t4) and 3 reads (r1,r2,r3).
# The corrispondance between read and table is: r1={t1=1,t2=1,t3=0,t4=1},r2={t1=0,t2=0,t3=1,t4=1},r3={t1=1,t2=0,t3=0,t4=1}
# The individual targets have prior = {t1=0.25,t2=0.25,t3=0.4,t4=0.1}
# We want to compute the posterior probability of each target given the reads r1,r2,r3 and the prior probability of each target.
# The posterior probability is computed as follows:
# P(t1|r1,r2,r3) = P(r1|t1) * P(r2|t1) * P(r3|t1) * P(t1) / P(r1) * P(r2) * P(r3)
# P(r1) * P(r2) * P(r3) is a constant and can be left out.
# Where (P(r1|t1) * P(r2|t1) * P(r3|t1)) cis the likelihood of the reads given the target t1.
# And P(t1) is the prior probability of the target t1.
# Same for the other targets.
# Computed the posterior probability of each target we can normalize them to 1 and use them as the new prior probability of each target.

# The computed posterior per target will be then used to set again the prior probability of each target and perform severla iterations of the EM algorithm.
# We start with a data frame of read_id and target_id pairs. Only possible pairs are shown therefore if a read_id does not map to a target it won't be represented in the data frame's rows.
# We divide the various stconvergence of the EM algorithm in different functions.

# Compute the intial prior as a function of single mappers reads per targets:
# Based on the single mappers it returns the prior per target:
# it takes in a matrix with single and multi mappers reads and returns dataframe with the prior per target
def compute_initial_prior(df):
    # Get the single mappers
    df_single = df[df['counts_per_read_id'] == 1]
    # Fix the index
    df_single = df_single.reset_index(drop=True)
    # Get the number of single mappers per target return a data frame with 2 columns target_id and count
    prior = df_single.groupby('target_id')['target_id'].count().reset_index(name='count_target')
    # Compute the prior as the number of single mappers per target divided by the total number of single mappers
    prior['prior'] = prior['count_target'] / prior['count_target'].sum()
    # Return only the target_id and the prior
    return prior[['target_id', 'prior']]

# Compute the posterior probability of each target given the whole set of multi-mappers and the current prior probability of each target:
# First compute the likelihood of the reads given the target tj.
# The it computes the posterior probability of each target given the reads r1,r2,r3 and the prior probability of each target.
def compute_posterior(df, current_prior):
    # Compute the likelihood:
    # Get the multi mappers
    df_multi = df[df['counts_per_read_id'] > 1]
    # Fix the index
    df_multi = df_multi.reset_index(drop=True)
    # Compute the likelihood as the inverse of the number of multi mappers per target
    df_multi.loc[:, 'likelihood'] = 1 / df_multi['counts_per_read_id']
    # Compute the likelihood per target
    likelihood = df_multi.groupby('target_id')['likelihood'].prod().reset_index(name='likelihood')
    # Compute the posterior:
    posterior = pd.merge(likelihood, current_prior, on='target_id')
    # Compute the posterior
    posterior['posterior'] = posterior['likelihood'] * posterior['prior']
    # Normalize the posterior
    posterior['posterior'] = posterior['posterior'] / posterior['posterior'].sum()
    # Return the posterior
    return posterior

# Perform the EM algorithm
# The function takes in input the dataframe with the multi mappers and single mapperS, the number of iterations, and the convergenceilon value.
# to make max_iter and convergence optional arguments we can set them to 1000 and 1e-6 respectively.
def EM(df,max_iter=1000,convergence=1e-6):
    # Compute the intial prior
    initial_prior = compute_initial_prior(df)
    current_prior = initial_prior.copy()
    # Start the iteration:
    # Diplay the progression as fraction of total iterations in the form of a progress bar
    for i in range(max_iter): # max_iter is the number of iterations
        # Display the progress
        sys.stderr.write('\rIteration {}/{}'.format(i + 1, max_iter))
        # Compute the posterior probability of each target given the whole set of multi-mappers and the current prior probability of each target
        posterior = compute_posterior(df, current_prior) 
        # Check the convergence of the algorithm
        if np.abs(posterior['posterior'] - current_prior['prior']).max() < convergence:
            break
        # Update the prior assign the posterior from posterior to the 'prior' column of the current_prior
        current_prior['prior'] = posterior['posterior'] 
    # Return the posterior
    return posterior


# Below is the main function
if __name__ == '__main__':
    # Parse the arguments
    args = parse_args()
    # Read the input file
    mat = read_input_file(args.input_file)
    # get unique raws based on read_id and target_id - as reads are unique entity they can map at most to one target
    # Later for ambiguous cases in which the same multimapper can map to the same target but in different positions those are going to be resolved
    mat_filter = mat.drop_duplicates(subset=['read_id', 'target_id'])
    # We compute the posterior probability of each target:
    posterior = EM(mat_filter, args.max_iter, args.convergence)
    # Normalize the posterior so that the final will sum to 1
    posterior['posterior'] = posterior['posterior'] / posterior['posterior'].sum()
    # use make_dir to create the output directory and save the posterior to the output directory
    make_dir(args.output_file)
    # Save the posterior save to the working directory the 'posterior_target_probabilities.txt' as tab delimited:
    posterior.to_csv(os.path.join(args.output_file, 'posterior_target_probabilities.txt'), sep='\t', index=False)
    # Starting from the original matrix we select raws with target_id in the posterior dataframe
    # We do this to then select the most probable target for each read. For cases where the same read maps to the same target in different positions we select randomly one of the positions and then the most probable given the posterior.
    # from the original mat we get all single mappers and all multi mappers
    mat_single = mat[mat['counts_per_read_id'] == 1]
    # we get the multi mappers from the original and select those that are in the posterior dataframe
    mat_multi = mat[mat['counts_per_read_id'] > 1]
    mat_multi = mat_multi[mat_multi['target_id'].isin(posterior['target_id'])]
    # Combine mat_multi and posterior
    mat_multi = pd.merge(mat_multi, posterior, on='target_id')
    # Shuffle the raws of the multi mappers - so that equal posterior per read_id are randomly selected
    mat_multi = mat_multi.sample(frac=1).reset_index(drop=True)
    # Sort by posterior so that higher posterior are first
    mat_multi = mat_multi.sort_values(by=['posterior'], ascending=False)
    # select the most probable target per read_id frop those with 0 posterior   
    mat_multi = mat_multi[mat_multi['posterior'] > 0].drop_duplicates(subset=['read_id'])
        
    #Combine single and multi mappers to obtain the final matrix with one target per read so all are going to be unique
    # Concatenate the two matrices selecting the first 13 columns
    mat_concat = pd.concat([mat_single.iloc[:, :3], mat_multi.iloc[:, :3]])
    # Sort the matrix by read_id
    mat_concat = mat_concat.sort_values(by=['read_id'])
    # Reset the index
    mat_concat = mat_concat.reset_index(drop=True)
    # Ensure one read_id is present only once
    mat_concat = mat_concat.drop_duplicates(subset=['read_id'])
    
    # Save the final matrix without column names
    mat_concat.to_csv(os.path.join(args.output_file, 'Final.bedpe'), sep='\t', index=False, header=False)
    