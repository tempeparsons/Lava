import math
import numpy as np

from scipy.stats import ttest_ind_from_stats

    
def remove_rows(df, grpsize1, grpsize2):

    cols = df.columns
 
    # fine
    # 4/4 1/4 
    # 3/4 1/4
    # 3/3 1/3 
    
    # borderline
    # 2/4 2/4
    
    # not fine
    # 2/4 1/4
    # 1/3 1/3 
    # 2/2 1/2 ?
    # 2/3 1/3 ?
    # 2/2 1/3 ?
    
    # Remove if up to half missing on both sides
    thresh1 =  grpsize1 // 2
    thresh2 =  grpsize2 // 2
   
    # Drop if less than half on both sides: 2/4 both sides is borderline and kept
    df = df.drop(df[(df.iloc[:,-1] <= thresh1) & (df.iloc[:,-2] < thresh2)].index)
    df = df.drop(df[(df.iloc[:,-1] < thresh1) & (df.iloc[:,-2] <= thresh2)].index)
    
    # Need at least 3 on one side to compare a singlular val
    df = df.drop(df[(df.iloc[:,-1] < 3) & (df.iloc[:,-2] == 1)].index)
    df = df.drop(df[(df.iloc[:,-1] == 1) & (df.iloc[:,-2] < 3)].index)
    
    # Any pure zeros are real zeros not NaN
    df.loc[df[cols[-2]] == 0, cols[:grpsize1]] = [0] * grpsize1
    df.loc[df[cols[-1]] == 0, cols[grpsize1:grpsize1+grpsize2]] = [0] * grpsize2
    
    # Any singlular dat has nobs set to 2 for t-test
    df.loc[df[cols[-2]] == 1, cols[-2]] = 2 # Nobs for t-test
    df.loc[df[cols[-1]] == 1, cols[-1]] = 2
    
    return df
    
    
def prpn_log2FCs_over_threshold(df, thresh_int):    
    fcs = df['grp1grp2_FC'].to_list()
    totalpos = len(([i for i in fcs if i > 0]))
    totalneg = len(([i for i in fcs if i < 0]))
    over_pos_thresh = len([i for i in fcs if i >= thresh_int])
    over_neg_thresh = len([i for i in fcs if i <= -thresh_int])
    prpn_pos = (over_pos_thresh/totalpos)*100 
    pos_pct = "{:.2f}".format(prpn_pos)
    prpn_neg = (over_neg_thresh/totalneg)*100
    neg_pct = "{:.2f}".format(prpn_neg)
    return pos_pct, neg_pct



def make_znorm(df):
    valuecols = df.columns.to_list()
    for col in valuecols:
        df[col] = (df[col] - df[col].mean())/df[col].std(ddof=0)
    df += 10
    return df



def ttest_from_stats_eqvar(df):

    tt_res = ttest_ind_from_stats(mean1 = df['zmean_grp1'], std1 = df['zstd_grp1'], nobs1 = df['znobs_grp1'], 
                                  mean2 = df['zmean_grp2'], std2 = df['zstd_grp2'], nobs2 = df['znobs_grp2'],
                                  equal_var=True)
    df['pvals'] = tt_res[1].tolist() #to_list
    df['Tpvals'] = -1*np.log2(tt_res[1].tolist()) #to_list

    tt_res = ttest_ind_from_stats(mean1 = df['zmean_grp1'], std1 = df['zstd_grp1_q95'], nobs1 = df['znobs_grp1_q95'], 
                                  mean2 = df['zmean_grp2'], std2 = df['zstd_grp2_q95'], nobs2 = df['znobs_grp2_q95'], 
                                  equal_var=True)
    df['pvals_q95'] = tt_res[1].tolist() #to_list
    df['Tpvals_q95'] = -1*np.log2(tt_res[1].tolist()) #to_list
    return df
    
    
def get_bins(pair_df):

    FCs = pair_df['grp1grp2_FC'].to_list()
    fcmin = np.nanmin(FCs)
    fcmax = np.nanmax(FCs)
    counts, edges = np.histogram(FCs, bins=50, range=(fcmin, fcmax))
    
    return(counts, edges)



def get_xlist(pair_df):
    counts, edges = get_bins(pair_df)

    xmin = None
    xmax = None

    for i in range(len(counts)-1):
        xmin = np.ceil(min(edges)-1)
        xmax = float(math.floor(max(edges)+1))

    xtemp = [abs(xmin), abs(xmax)]
    biggest = max(xtemp)
    xlist = [-biggest, biggest]
    
    return xlist



def get_splitxlist(pair_df):
    counts, edges = get_bins(pair_df)
    xmiddles = []

    x0 = None
    x1 = None
    x2 = None
    x3 = None
    x4 = None
    x5 = None

    xmin = None
    xmax = None

    for i in range(len(counts)-1):
        if counts[i] != 0 and counts[i+1] == 0 and counts[i+2] == 0 and counts[i+3] == 0:
            xmiddles.append(edges[i+3])
        else:
            xmin = np.ceil(min(edges)-1)
            xmax = float(math.floor(max(edges)+1))

    x0 = np.ceil(min(edges)-1)
    x1 = np.ceil(xmiddles[0]-1)
    x2 = np.ceil(-xmiddles[1])
    x3 = float(math.floor(xmiddles[1]))
    x4 = float(math.floor(abs(xmiddles[0])+1))
    x5 = float(math.floor(max(edges)+1))

    splitxlist = [x0, x1, x2, x3, x4, x5]
    splitxlist = sorted(splitxlist)
    
    return splitxlist
    
    
