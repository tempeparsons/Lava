
from matplotlib import pyplot as plt

def plot_sums_means_nans(df, value_cols):
    totals = df[value_cols].sum(axis=0)
    means = df[value_cols].mean(axis=0)
    nancounts = df[value_cols].isna().sum()
    contents = [totals, means, nancounts]
    titles = ['column totals', 'column means', 'column nan counts']
    
    fig, axs = plt.subplots(3, figsize=(6,3), sharex=True)
    
    for i, ax in enumerate(axs.flat):
        ax.bar(value_cols, contents[i])
        ax.set_title(titles[i], fontsize=10)
        ax.margins(x=0)
        ax.tick_params(axis='x', rotation=90)
        
    plt.subplots_adjust(hspace=0.4)
    
    # # # # 
    
    plt.show()
    
    
def boxplots(df, value_cols):
    fig, axs = plt.subplots(1,1, figsize=(8,4))
    df.boxplot(column=value_cols, grid=False, rot=45, fontsize=8)     
    
    
    
def histograms(df, value_cols):
    #cols = df.columns[value_cols].to_list()
    params = {'axes.titlesize':'8',
              'xtick.labelsize':'7',
              'ytick.labelsize':'7', 
              'axes.titley': '0.7'}
    plt.rcParams.update(params)
    df.hist(column=value_cols, 
            alpha=0.5, 
            bins=50, 
            grid=True, 
            sharex=True, sharey=True, 
            xlabelsize=6, ylabelsize=6,
            figsize=(9,9))
    plt.subplots_adjust(hspace=0.1, wspace=0.1) 
    plt.show()

    
    
def xy_plots(df, value_cols, ncols):
    nrows = df.shape[1] // ncols + (df.shape[1] % ncols > 0)
    cols = value_cols
    fig, axs = plt.subplots(nrows,ncols, figsize=(8,8), sharex=True, sharey=True)
    for i in range(len(cols)-1):
        axs.flat[i].scatter(df[cols[i]], df[cols[i+1]], s=5, alpha=0.6)
        axs.flat[i].set_title(f'{cols[i]} vs {cols[i+1]}', fontsize=8)
        ymin, ymax = axs.flat[i].get_ylim()
        axs.flat[i].axline([ymax*0.5,ymax*0.5], slope=1, c='r', lw=0.5)    
    plt.subplots_adjust(wspace=0.05, hspace=0.2)


    
def remove_rows(df, grpsize1, grpsize2):
    cols = df.columns
    if grpsize1 == 4 and grpsize2 == 4:
        df = df.drop(df[(df.iloc[:,-1] == 2) & (df.iloc[:,-2] == 1)].index)
        df = df.drop(df[(df.iloc[:,-1] == 1) & (df.iloc[:,-2] == 2)].index)
        df = df.drop(df[(df.iloc[:,-1] == 0) & (df.iloc[:,-2] == 2)].index)
        df = df.drop(df[(df.iloc[:,-1] == 2) & (df.iloc[:,-2] == 0)].index)
        df = df.drop(df[(df.iloc[:,-1] == 1) & (df.iloc[:,-2] == 1)].index)
        df = df.drop(df[(df.iloc[:,-1] == 0) & (df.iloc[:,-2] == 1)].index)
        df = df.drop(df[(df.iloc[:,-1] == 1) & (df.iloc[:,-2] == 0)].index)
        df = df.drop(df[(df.iloc[:,-1] == 0) & (df.iloc[:,-2] == 0)].index)
        df.loc[df[cols[-2]] == 0, cols[:4]] = [0,0,0,0]
        df.loc[df[cols[-1]] == 0, cols[4:8]] = [0,0,0,0]
        df.loc[df[cols[-2]] == 1, cols[-2]] = 2
        df.loc[df[cols[-1]] == 1, cols[-1]] = 2
        return df
    if grpsize1 == 3 or grpsize2 == 3:
        df = df.drop(df[(df.iloc[:,-1] == 1) & (df.iloc[:,-2] == 1)].index)
        df = df.drop(df[(df.iloc[:,-1] == 0) & (df.iloc[:,-2] == 1)].index)
        df = df.drop(df[(df.iloc[:,-1] == 1) & (df.iloc[:,-2] == 0)].index)
        df = df.drop(df[(df.iloc[:,-1] == 0) & (df.iloc[:,-2] == 0)].index)
        if grpsize2 == 3:
            df.loc[df[cols[-2]] == 0, cols[:4]] = [0,0,0,0]#######################
            df.loc[df[cols[-1]] == 0, cols[4:7]] = [0,0,0]
        if grpsize1 == 3:
            df.loc[df[cols[-2]] == 0, cols[:3]] = [0,0,0]
            df.loc[df[cols[-1]] == 0, cols[3:7]] = [0,0,0,0]
        if grpsize1 == 3 and grpsize2 == 4:
            df.loc[df[cols[-4]] == 0, cols[:3]] = [0,0,0]
            df.loc[df[cols[-1]] == 0, cols[3:6]] = [0,0,0]
        df.loc[df[cols[-2]] == 1, cols[-2]] = 2
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



def pvalFC_hists(plotdictkey):
    df = plotdict[plotdictkey]
    fig, axs = plt.subplots(1,2, figsize = (8,2), sharex='col', sharey=False)
    pvals = df['pvals'].to_list()
    FCs = df['grp1grp2_FC'].to_list()
    axs[0].hist(pvals, bins=50)
    axs[0].set_title(f'p-value distrbution for {plotdictkey}', fontsize=6)
    axs[1].hist(FCs, bins=50)
    axs[1].set_title(f'fold change distribution for {plotdictkey}', fontsize=6)
    axs[0].set_ylabel('frequency', fontsize=6)
    axs[0].set_xlabel('value', fontsize=6)
    axs[1].set_xlabel('value', fontsize=6)

    
    
    
def get_bins(plotdictkey):
    df = plotdict[plotdictkey]
    FCs = df['grp1grp2_FC'].to_list()
    fcmin = np.nanmin(FCs)
    fcmax = np.nanmax(FCs)
    counts, edges = np.histogram(FCs, bins=50, range=(fcmin, fcmax))
    return(counts, edges)



def get_xlist(df_pair):
    counts, edges = get_bins(df_pair)

    xmin = None
    xmax = None

    for i in range(len(counts)-1):
        xmin = np.ceil(min(edges)-1)
        xmax = float(math.floor(max(edges)+1))

    xtemp = [abs(xmin), abs(xmax)]
    biggest = max(xtemp)
    xlist = [-biggest, biggest]
    
    return xlist



def get_splitxlist(df_pair):
    counts, edges = get_bins(df_pair)
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
