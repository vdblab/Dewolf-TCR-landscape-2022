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
    args = parser.parse_args()
    return args

def filter_tcrs():
  # remove all out of frame one
  pass
def norm_jsdiv(P, Q):
  pass

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

def pivot_tcr(both, bythis="amino_acid"):
  thisagg = "sum" 
  return (both[["sample_name", bythis,"productive_frequency"]]
    .pivot_table(
      index="sample_name", 
      columns=bythis,
      values="productive_frequency",
      aggfunc=thisagg,
      fill_value=0)
    )
  
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
  Adf = pd.read_csv(args.A, sep="\t", usecols=cols_of_interest).query('frame_type == "In"')
  Bdf = pd.read_csv(args.B, sep="\t", usecols=cols_of_interest).query('frame_type == "In"')
  both = pd.concat([Adf, Bdf], axis=0)
  Asize = Adf["productive_templates"][0]
  Bsize = Bdf["productive_templates"][0]
  print(Asize, Bsize)
  depth_threshold = min(Asize, Bsize)

  both_filt = (both
    .assign(lowdepth_templates = lambda x: (x.templates * depth_threshold) / x.productive_templates)
    .query('lowdepth_templates >= 1')
  )
  if verbose: 
    print(both_filt.head())
    print(both.shape)
    print(both_filt.shape)
  new_templates_totals = both_filt[["sample_name", "lowdepth_templates"]].groupby("sample_name").sum())

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

  
  print(depth_threshold)
  both_w_aa = pivot_tcr(both, bythis="amino_acid")
  both_w_nt = pivot_tcr(both, bythis="rearrangement")
  both_w_aa_filt = both_w_aa
  both_w_nt_filt = both_w_nt

  # this is way faster than manually calculating (eg .0015 vs .052s for the balb_1_blood vs balb_2_blood)
  jsd_nt_raw = distance.jensenshannon(both_w_nt.iloc[0], both_w_nt.iloc[1], base=2) ** 2
  # NT JSD Filter
  jsd_nt_filt = distance.jensenshannon(both_w_nt.iloc[0], both_w_nt.iloc[1], base=2) ** 2
  # AA JSD raw
  jsd_aa_raw = distance.jensenshannon(both_w_aa.iloc[0], both_w_aa.iloc[1], base=2) ** 2
  # AA JSD Filter
  jsd_aa_filt = distance.jensenshannon(both_w_aa.iloc[0], both_w_aa.iloc[1], base=2) ** 2
  print(jsd_nt_raw, jsd_nt_filt, jsd_aa_raw, jsd_aa_filt)
  
if __name__ == "__main__":
  args=get_args()
  main(args)
