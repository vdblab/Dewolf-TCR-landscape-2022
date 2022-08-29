#!/usr/bin/env python3
from scipy.spatial import distance
# 
# 
# 
import numpy as np
import pandas as pd
import argparse
import os
import time
import sys

__version__ = "0.0.1"
def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="compute JSD and normalized JSD on TCR sequencing files from AdaptiveBiotech " +
        "version " + __version__,
        add_help=False)
    parser.add_argument("A", action="store",
                        help="First TCR seq file")
    parser.add_argument("B", action="store",
                        help="Second TCR seq file")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="print debugging info")
    parser.add_argument("-d", "--header", action="store_true",
                        help="print header info")
    args = parser.parse_args()
    return args

def filter_tcrs():
  # remove all out of frame one
  pass
def norm_jsdiv(P, Q):
  pass


def morisita(x, y):
  # taken verbatim from abdiv
    Nx = sum(x)
    Ny = sum(y)
    lambda_x = sum(x**2)/(Nx**2)
    lambda_y = sum(y**2)/(Ny**2)
    m = 1 - 2 * sum(x * y)/((lambda_x + lambda_y) * Nx * Ny)
    return(m)

def jsdiv(P, Q):
  
    """Compute the Jensen-Shannon divergence between two probability distributions.
    # https://stackoverflow.com/questions/15880133/jensen-shannon-divergence
    https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
    Input
    -----
    P, Q : array-like
        Probability distributions of equal length that sum to 1
    """

    def _kldiv(A, B):
        return np.sum([v for v in A * np.log2(A/B) if not np.isnan(v)])

    P = np.array(P)
    Q = np.array(Q)

    M = 0.5 * (P + Q)

    return 0.5 * (_kldiv(P, M) +_kldiv(Q, M))
  
# depth_threshold <- min(tcr_corrs_equalized$productive_templates.x[i], tcr_corrs_equalized$productive_templates.y[i] )
#   tcr_depth_filtered <- tcr_raw %>% 
#     filter(sample_name %in% these_samples) %>% 
#     select(sample_name, rearrangement, amino_acid, templates, productive_templates, productive_frequency) %>%
#     mutate(lowdepth_templates = (templates *depth_threshold) / productive_templates) %>%
#     filter(lowdepth_templates >=1) %>%
#     group_by(sample_name) %>%
#     mutate(total_lowdepth_templates = sum(lowdepth_templates)) %>%
#     ungroup() %>%
#     mutate(normalized_lowdepth_templates = (lowdepth_templates * min(total_lowdepth_templates)) / total_lowdepth_templates) %>%
#     mutate(lowdepth_productive_frequency = normalized_lowdepth_templates / min(total_lowdepth_templates)) %>%
#     data.table::as.data.table() %>% 
#     data.table::dcast(sample_name~rearrangement, fill = 0, value.var="lowdepth_productive_frequency" ) %>%
#     column_to_rownames("sample_name") %>% 
#     data.matrix()

def pivot_tcr(both, bythis="amino_acid", fill=True):
  # using speed up from https://stackoverflow.com/questions/55404617/faster-alternatives-to-pandas-pivot-table
  thisagg = "sum" 
  # return (both[["sample_name", bythis,"productive_frequency"]]
  #   .pivot_table(
  #     index="sample_name", 
  #     columns=bythis,
  #     values="productive_frequency",
  #     aggfunc=thisagg,
  #     fill_value=0 if fill else np.nan)
  #   )
  widedf = both[["sample_name", bythis,"productive_frequency"]].groupby(["sample_name", bythis]).agg({"productive_frequency":thisagg}).unstack(level=bythis)
  if fill:
    return(widedf.fillna(0))
  else:
    return(widedf)
  
def main(args):
  verbose = args.verbose
  cols_of_interest = [
    "sample_name",
    "frame_type",
    "amino_acid",
    "rearrangement",
    "productive_frequency",
    "templates", 
    "productive_templates"
  ]
  # still useful to calculate some things on same input file pairs (eg "overlapping pairs"")
  same_file = False
  if (args.A == args.B):
    same_file = True
  Adf = pd.read_csv(args.A, sep="\t", usecols=cols_of_interest).query('frame_type == "In"')
  Bdf = pd.read_csv(args.B, sep="\t", usecols=cols_of_interest).query('frame_type == "In"')
  both = pd.concat([Adf, Bdf], axis=0)
  Asize = Adf["productive_templates"].iloc[0]
  Bsize = Bdf["productive_templates"].iloc[0]
  if verbose: print(Asize, Bsize)
  depth_threshold = min(Asize, Bsize)

  both_filt = (both
    .assign(lowdepth_templates = lambda x: (x.templates * depth_threshold) / x.productive_templates)
    .query('lowdepth_templates >= 1')
  )
  if verbose: 
    print(both_filt.head())
    print(both.shape)
    print(both_filt.shape)
  new_templates_totals = both_filt[["sample_name", "lowdepth_templates"]].groupby("sample_name").sum()
    # filter(lowdepth_templates >=1) %>%
    # group_by(sample_name) %>%
    # mutate(total_lowdepth_templates = sum(lowdepth_templates)) %>%
    # ungroup() %>%
    # mutate(normalized_lowdepth_templates = (lowdepth_templates * min(total_lowdepth_templates)) / total_lowdepth_templates) %>%
    # mutate(lowdepth_productive_frequency = normalized_lowdepth_templates / min(total_lowdepth_templates)) %>%
    # data.table::as.data.table() %>%
    # data.table::dcast(sample_name~rearrangement, fill = 0, value.var="lowdepth_productive_frequency" ) %>%
    # column_to_rownames("sample_name") %>%
    # data.matrix()

  
  if verbose: print(depth_threshold)
  if not same_file:
    both_w_aa = pivot_tcr(both, bythis="amino_acid")
    both_w_nt = pivot_tcr(both, bythis="rearrangement")
    both_w_aa_filt = pivot_tcr(both_filt, bythis="amino_acid" )
    both_w_nt_filt = pivot_tcr(both_filt, bythis="rearrangement")
    # this is way faster than manually calculating (eg .0015 vs .052s for the balb_1_blood vs balb_2_blood)
    jsd_nt_raw = distance.jensenshannon(both_w_nt.iloc[0], both_w_nt.iloc[1], base=2) ** 2
    # AA JSD raw
    jsd_aa_raw = distance.jensenshannon(both_w_aa.iloc[0], both_w_aa.iloc[1], base=2) ** 2
    # NT Morisita raw 
    morisita_nt_raw = morisita(x=both_w_nt.iloc[0], y=both_w_nt.iloc[1])
  
    # this can occur when filtering pairs with extreme differences in size
    # filtering/normalizing removes all tempaltes from the more diverse sample
    if both_w_nt_filt.shape[0] !=2:
      jsd_nt_filt = np.nan
      jsd_aa_filt = np.nan
      new_templates_totals_tup = (np.nan, np.nan)
    else :
      jsd_nt_filt = distance.jensenshannon(both_w_nt_filt.iloc[0], both_w_nt_filt.iloc[1], base=2) ** 2
      jsd_aa_filt = distance.jensenshannon(both_w_aa_filt.iloc[0], both_w_aa_filt.iloc[1], base=2) ** 2
      new_templates_totals_tup = ( new_templates_totals['lowdepth_templates'][0] ,  new_templates_totals['lowdepth_templates'][1] )
    # overlapping clones 
    overlapping_clones = pivot_tcr(both, bythis="rearrangement", fill=False).dropna(axis="columns", how='any').shape[1]
    overlapping_aa = pivot_tcr(both, bythis="amino_acid", fill=False).dropna(axis="columns", how='any').shape[1]
  
  else:
    jsd_nt_raw = 0
    jsd_aa_raw = 0
    jsd_nt_filt = 0
    jsd_aa_filt = 0
    morisita_nt_raw = 0
    new_templates_totals_tup = (np.nan, np.nan)
    overlapping_clones = Adf["rearrangement"].nunique()
    overlapping_aa  = Adf["amino_acid"].nunique()
  if args.header:
    print(f"fileA\tfileB\tsizeA\tsizeB\tnorm_sizeA\tnorm_sizeB\tjsd_nt_raw\tjsd_nt_norm\tjsd_aa_raw\tjsd_aa_norm\tmorisita_nt_raw\toverlapping_clones\t{overlapping_aa}")
    
  print(f"{os.path.basename(args.A).replace('.tsv', '')}\t{os.path.basename(args.B).replace('.tsv', '')}\t{Asize}\t{Bsize}\t{new_templates_totals_tup[0]}\t{new_templates_totals_tup[1]}\t{jsd_nt_raw}\t{jsd_nt_filt}\t{jsd_aa_raw}\t{jsd_aa_filt}\t{morisita_nt_raw}\t{overlapping_clones}\t{overlapping_aa}")

if __name__ == "__main__":
  args=get_args()
  try:
    main(args)
  except Exception as e:
    sys.stderr.write(", ".join(sys.argv))
    raise(e)
  
##  parallel AB-correlations.py {1} {2} > tmp :::: manifest.txt ::::+ manifest.txt
    
  
