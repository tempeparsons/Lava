import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patheffects as PathEffects

DPI = 300


def plot_sums_means_nans(df, value_cols, pdf, colors=['#D00000','#A0A000','#0080FF'], edge_color='#404040'):
    
    totals = df[value_cols].sum(axis=0)
    means = df[value_cols].mean(axis=0)
    nancounts = df[value_cols].isna().sum()
    contents = [totals, means, nancounts]
    n_cols = len(value_cols)
    
    titles = ['Totals', 'Averages', 'Missing values']
    y_labels = ['Column total', 'Column mean', 'Column NaN count']
    
    fig, axs = plt.subplots(3, figsize=(8,8), sharex=False)
    x_pad = 0.6
    x_vals = np.arange(0, n_cols)
    
    for i, ax in enumerate(axs.flat):
        ax.bar(x_vals, contents[i], color=colors[i], edgecolor=edge_color)
        ax.set_title(titles[i], fontsize=12, color=colors[i], fontdict={'fontweight':'bold'} )
        ax.set_ylabel(y_labels[i])
        ax.set_xlim(-x_pad, n_cols-1+x_pad)
        ax.set_xticks(x_vals)
        ax.set_xticklabels(1+x_vals)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        for j, col in enumerate(value_cols):
            txt = ax.text(j, 0.1*max(contents[i]), col, va='bottom', ha='center', rotation=90, color='#FFFFFF')
            txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground=edge_color)])

    ax.set_xlabel('Data column/sample')    
    plt.subplots_adjust(hspace=0.35, bottom=0.1, top=0.95, right=0.95)
    
    if pdf:
      pdf.savefig(dpi=DPI)
    else:
      plt.show()
    
    plt.close()
    
    
def boxplots(df, value_cols, pdf, scatter_colors=['#D00000','#A0A000','#0080FF'], scatter_width=0.7):
    
    box_style = {'linewidth':1, 'color':'#0080FF', 'markerfacecolor':'#D00000'}
    line_style = {'linewidth':1, 'color':'#0080FF'}
    outlier_style = {'linewidth':1, 'markeredgecolor':'#D00000', 'markerfacecolor':'#FF8080'}
    med_style = {'linewidth':1, 'color':'#000000'}
    
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,8))
    df.boxplot(column=value_cols, ax=ax1, grid=False, rot=45, fontsize=8, boxprops=box_style,
               whiskerprops=line_style, capprops=line_style, flierprops=outlier_style, medianprops=med_style)     
    ax1.set_title('Sample value distributions')

    ns = len(scatter_colors)
    nc = len(value_cols)
    
    for i, col in enumerate(value_cols):
      y_vals = np.sort(df[col])
      
      # Random but avoidung clumps
      dx = np.random.uniform(scatter_width/4, scatter_width/2, len(y_vals))
      x_vals = i + (np.cumsum(dx) % scatter_width) - scatter_width/2
      
      color = scatter_colors[i % ns]
      ax2.scatter(x_vals, y_vals, color=color, alpha=0.4, s=5)
    
    x_ticks = np.arange(0, nc)   
    ax2.set_xlim(-0.5, nc-0.5)  
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_ticks+1)
    ax2.set_xlabel('Sample column')
    ax2.tick_params('x', top=True)
    
    plt.subplots_adjust(hspace=0.35, top=0.90, bottom=0.1, left=0.05, right=0.95)    
    
    if pdf:
      pdf.savefig(dpi=DPI)
    else:
      plt.show()
    
    plt.close()
    
    
def histograms(df, value_cols, pdf):
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
    
    plt.subplots_adjust(hspace=0.1, wspace=0.1, top=0.95, bottom=0.05, left=0.05, right=0.95) 
    
    if pdf:
      pdf.savefig(dpi=DPI)
    else:
      plt.show()
    
    plt.close()
    

    
    
def xy_plots(df, value_cols, ncols, pdf):
    nrows = df.shape[1] // ncols + (df.shape[1] % ncols > 0)
    cols = value_cols
    fig, axs = plt.subplots(nrows, ncols, figsize=(8,8), sharex=True, sharey=True)
    for i in range(len(cols)-1):
        axs.flat[i].scatter(df[cols[i]], df[cols[i+1]], s=5, alpha=0.6)
        txt = axs.flat[i].set_title(f'{cols[i]} vs {cols[i+1]}', fontsize=8)
        txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
        ymin, ymax = axs.flat[i].get_ylim()
        axs.flat[i].axline([ymax*0.5,ymax*0.5], slope=1, c='r', lw=0.5)    
    
    # switch off unused axes
    for i in range(len(value_cols)-1, nrows*ncols):
        axs.flat[i].set_visible(False)
        
    plt.subplots_adjust(wspace=0.05, hspace=0.2, top=0.95, bottom=0.05, left=0.05, right=0.95)
    
    if pdf:
      pdf.savefig(dpi=DPI)
    else:
      plt.show()
    
    plt.close()

    
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
