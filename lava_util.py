import math
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patheffects as PathEffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import ttest_ind_from_stats

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
        ax.bar(x_vals, contents[i], color=colors[i], edgecolor=edge_color, alpha=0.75)
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
    
    #params = {'axes.titlesize':'8',
    #          'xtick.labelsize':'7',
    #          'ytick.labelsize':'7', 
    #          'axes.titley': '0.7'}
              
    #plt.rcParams.update(params)
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



def pvalFC_hists(plotdict, pdf, fontsize=8, colors=['#D00000','#0080FF'], nbins=50):
    
    ###  Splt this between comparison sets
    
    keys = sorted(plotdict.keys())
    n_rows = len(keys)
    
    fig, axs = plt.subplots(n_rows, 2, figsize = (8, n_rows*1.2), sharex='col')
    
    c0, c1 = colors
    x0_min = -1.0/nbins
    x0_max = 1.0 - x0_min
     
    for row, key in enumerate(keys):
         df = plotdict[key]
         ax0 = axs[row, 0]
         ax1 = axs[row, 1]
         
         pvals = df['pvals'].to_list()
         FCs = df['grp1grp2_FC'].to_list()
         
         if row == 0:
             ax0.set_title(f'P-value distrbutions', color=c0)    
             ax1.set_title(f'Fold change distributions', color=c1)
 
         txt0 = ax0.text(0.05, 0.8, key, transform=ax0.transAxes, fontsize=12)
         txt1 = ax1.text(0.05, 0.8, key, transform=ax1.transAxes, fontsize=12)
         txt0.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
         txt1.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
         
         #ax0.hist(pvals, bins=nbins, color=c0, range=(0.0, 1.0))
         hist, edges = np.histogram(pvals, bins=nbins, range=(0.0, 1.0))
         x_vals = 0.5 * (edges[:-1] + edges[1:])
         ax0.plot(x_vals, hist, color=c0, linewidth=1)
         ax0.fill_between(x_vals, hist, 0, color=c0, alpha=0.5)

         ax0.set_xlim(x0_min, x0_max)
          
         #ax1.hist(FCs, bins=nbins, color=c1)
         hist, edges = np.histogram(FCs, bins=nbins)
         x_vals = 0.5 * (edges[:-1] + edges[1:])
         ax1.axvline(0.0, linewidth=1, linestyle='--', color='#808080', alpha=0.5)
         ax1.plot(x_vals, hist, color=c1, linewidth=1)
         ax1.fill_between(x_vals, hist, 0, color=c1, alpha=0.5)
         
         ax0.set_ylabel('Bin count', fontsize=fontsize)
         
         ax0.set_ylim(0)
         ax1.set_ylim(0)
         
         if row == (n_rows-1):
             ax0.set_xlabel('p-value', fontsize=fontsize)
             ax1.set_xlabel('$log_2$ fold-change', fontsize=fontsize)
         #else:
         #    ax0.set_xticklabels([])
         #    ax1.set_xticklabels([])
         
    
    plt.subplots_adjust(wspace=0.15, hspace=0.1, top=0.95, bottom=0.05, left=0.1, right=0.95)

    if pdf:
      pdf.savefig(dpi=DPI)
    else:
      plt.show()
    
    plt.close()
    
    
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
    
    
def volcano(pdf, pair_name, df, q95, FClim, pthresh, colors, hq_only=False, hit_labels=True, print_hits=True, markers=None, lw=0.25, ls='--'):
    
    group1, group2 = pair_name.split(':')
    q951, q952 = q95
    
    if not markers:
      markers = []
      
    datacols = [df.index, df['grp1grp2_FC'], df['Tpvals'], df['Tpvals_q95'], df['zstd_grp1_q95'], df['zstd_grp2_q95']]

    if hq_only:
        plotlist = datacols[:3]
    else:
        plotlist = datacols[:6]
    
    pos_overFC = []
    neg_overFC = []
    
    color_pos, color_neg, color_low = colors
    xlims = get_xlist(df)
 
    fig, ax = plt.subplots(figsize = (10,8), sharey=True)
 
    ax.set_xlim(xlims[0], xlims[1])

    ax.axvline(FClim,  ls=ls, c='k', alpha=0.5, lw=lw)
    ax.axvline(-FClim, ls=ls, c='k', alpha=0.5, lw=lw)
    ax.axhline(pthresh,ls=ls, c='k', alpha=0.5, lw=lw)

    ax.set_title(pair_name, fontsize=18)
    ax.set_ylabel('-log2 transformed p-value')
    ax.set_xlabel('log2 fold change')
 
    ax.scatter(plotlist[1], plotlist[2], s=8, color='w', edgecolors='darkgrey', lw=lw)

    marker_text = []

    for i, name in enumerate(plotlist[0]):

        if name in markers:
            marker_text.append((plotlist[1][i], plotlist[2][i], name))
            ax.scatter(plotlist[1][i], plotlist[2][i], s=3, color='k')

        if plotlist[1][i] >= FClim and plotlist[2][i] >= pthresh:
            pos_overFC.append(name)
            ax.scatter(plotlist[1][i], plotlist[2][i], s=10, color=color_pos)
            if hit_labels:
                ax.annotate(name, xy=(plotlist[1][i], plotlist[2][i]), fontsize=6)

        if plotlist[1][i] <= -FClim  and plotlist[2][i] >= pthresh:
            neg_overFC.append(name)
            ax.scatter(plotlist[1][i], plotlist[2][i], s=10, color=color_neg)
            if hit_labels:
                ax.annotate(name, xy=(plotlist[1][i], plotlist[2][i]), fontsize=6)

        if not hq_only:
            if plotlist[3][i] == q951 or plotlist[4][i] == q952:
                ax.scatter(plotlist[1][i], plotlist[2][i], s=6, color=color_low)

    for tx in marker_text:
            ax.annotate(tx[2], xy=(tx[0], tx[1]),  xytext=(tx[0], tx[1]),
                    arrowprops=dict(arrowstyle='-', fc="k", ec="k", lw=0.5, relpos=(0.25, 0.5)),
                    bbox=dict(pad=-2, facecolor="none", edgecolor="none"),
                    ha="left", va="center", size=6)
 
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=f'Greater in {group1}',markerfacecolor=color_pos, markersize=8, linestyle=''),
                       Line2D([0], [0], marker='o', color='w', label=f'Greater in {group2}',markerfacecolor=color_neg, markersize=8, linestyle='')]
    ax.legend(handles=legend_elements)
    
    if print_hits:
        print('greater in grp1:', pos_overFC)
        print('greater in grp2:', neg_overFC)

    if pdf:
      pdf.savefig(dpi=DPI)
    else:
      plt.show()
    
    plt.close()
    
