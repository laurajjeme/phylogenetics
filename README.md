# phylogenetics

## tip2tip.pl

Will calculates for each taxon the average tip-to-tip distance between this given taxon and the 50% 'most distant' other taxa. The reason behind this is to avoid a bias for parts of the trees that more heavily sampled. It then looks for 'outliers' in the distribution of these values.

The definition for this is calculated as follows:
   a) Find the value of the 25th percentile (x) and the 75th percentile (y) of the average tip-to-tip  branch-length distribution
   b) Calculate the difference fs = y - x
   c) Looks for any tip-to-tip branch length that is over the following values:
       if length > y + (3 * fs)     ==> will be flagged as an extreme outlier
       if length > y + (1.5 * fs)   ==> will be flagged as an outlier
       
Made in collaboration with Dr. Ed Susko.
