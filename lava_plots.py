import math
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.patheffects as PathEffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.cluster import hierarchy
from scipy.spatial import distance

import lava_util as util

DPI = 300

def _watermark(fig, text=None):
    
    if not text:
      text = f'Lava v{VERSION}'
      
    pad = 0.1/fig.get_figwidth()
    
    plt.text(1.0-pad, pad, text, transform=fig.transFigure, color='#000000', alpha=0.25, va='bottom', ha='right', fontsize=8.0)

    

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
    _watermark(fig)
    
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
    _watermark(fig)
    
    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()
    
    
def histograms(df, value_cols, pdf, colors=['#D00000','#A0A000','#0080FF'], nbins=50):

    cols = list(df[value_cols])
    n_cols = len(cols)
    
    nx = int(math.ceil(math.sqrt(n_cols)))
    ny = int(math.ceil(n_cols/nx))
    
    fig, axarr = plt.subplots(ny, nx, squeeze=False, figsize=(9,9), sharex=True, sharey=True)
     
    for i in range(nx*ny):
        ax = axarr.flat[i]
 
        if i < n_cols:
            color = colors[i % len(colors)]
            data = df[value_cols[i]]
            data = data[~np.isnan(data)]
            hist, edges = np.histogram(data[~np.isnan(data)], bins=nbins)
            
            x_vals =  0.5 * (edges[:-1] + edges[1:])
            ax.plot(x_vals, hist, color=color, linewidth=1.0)
            ax.fill_between(x_vals, hist, 0, color=color, alpha=0.5) 
            txt = ax.text(0.05, 0.9, cols[i], transform=ax.transAxes, fontsize=12)
            txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
            
            n = len(data)
            #mu = np.mean(data)
            #sig = np.std(data)
            t = f'n={n:,}' #'\n$\mu$={mu:.2}\n$\sigma$={sig:.2}'
            
            txt = ax.text(0.95, 0.9, t, ha='right', va='top', color=color, transform=ax.transAxes, fontsize=8)
            txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
            
            if (i // nx) == ny-1:
                ax.set_xlabel('Data value')
            
            if (i % nx) == 0:
                ax.set_ylabel('Bin count')
            
        else:
            ax.set_visible(False)

    ax.set_ylim(0.0)
    
    plt.subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.1, left=0.1, right=0.95) 
    _watermark(fig)
    
    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()
        
    
def xy_plots(df, value_cols, ncols, pdf, colors=['#0080FF','#A0A000','#D00000'], figsize=10.0, bg_color='#B0B0B0'):
    
    cmap2 = LinearSegmentedColormap.from_list('gud', [bg_color, colors[0]])
    cmap1 = LinearSegmentedColormap.from_list('med', [bg_color, colors[1]])
    cmap0 = LinearSegmentedColormap.from_list('bad', [bg_color, colors[2]])
    
    cols = list(df[value_cols])
    n = len(cols)    
    
    fig, axarr = plt.subplots(n, n, squeeze=False, figsize=(figsize,figsize), sharex=True, sharey=True)
    
    gridsize = int(400//n)
    
    fontsize = min(1e3/(figsize * n), 14)
     
    for i in range(n):
        data1 = df[value_cols[i]]
        
        for j in range(n):
            ax = axarr[i,j]
            k = (i*n+j)

            if i != j:
                data2 = df[value_cols[j]]
                valid = ~(np.isnan(data1) | np.isnan(data2))
                color = colors[k % len(colors)]
                rho = np.corrcoef([data1[valid], data2[valid]])[0,1]
                
                if rho > 0.9:
                    cmap = cmap2
                elif rho > 0.7:
                    cmap = cmap1
                else:
                    cmap = cmap0
                
                ax.hexbin(data1[valid], data2[valid], cmap=cmap, gridsize=gridsize, mincnt=1, edgecolors='none')

                if i > j:
                    txt = ax.text(0.5, 0.5, f'{rho:.2}', ha='center', va='center', transform=ax.transAxes, fontsize=1.5*fontsize)
                    txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
 
            ax.set_xticks([])
            ax.set_yticks([])
 
            if i == (n-1):
                ax.set_xlabel(f'{cols[j]}', fontsize=fontsize)
                
            elif i == 0:    
                ax.xaxis.set_label_position('top') 
                ax.set_xlabel(f'{cols[j]}', fontsize=fontsize)
  
            if j == 0:
                ax.set_ylabel(f'{cols[i]}', fontsize=fontsize)
                
            elif j == (n-1):
                ax.yaxis.set_label_position('right')     
                ax.set_ylabel(f'{cols[i]}', fontsize=fontsize)

    plt.subplots_adjust(wspace=0.0, hspace=0.0, top=0.95, bottom=0.05, left=0.05, right=0.95)
    _watermark(fig)
    
    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()


def correlation_plot(df, value_cols, pdf, colors=['#0080FF','#A0A000','#D00000'], figsize=10.0, bg_color='#B0B0B0'):
    
    cmap = LinearSegmentedColormap.from_list('gud', colors)
    
    cols = list(df[value_cols])
      
    fig = plt.figure()
    fig.set_size_inches(figsize,figsize)
 
    # Matrix
    ax0 = fig.add_axes([0.15, 0.15, 0.70, 0.70]) # x, y, w, h
    # Right Dendrogram
    ax1 = fig.add_axes([0.85, 0.15, 0.10, 0.70])
    # Top Dendrogram
    ax2 = fig.add_axes([0.15, 0.85, 0.70, 0.10])
    ax2.set_title('Sample correlation matrix')
    # Colorbar
    ax3 = fig.add_axes([0.87, 0.85, 0.05, 0.10])
    
    data = np.nan_to_num(np.array([df[col] for col in cols]))
    
    corr_mat = distance.pdist(data, metric='correlation') # Correlation 'distance" : 1.0 - rho
    
    linkage = hierarchy.ward(corr_mat)
    order = hierarchy.leaves_list(linkage)[::-1]
    
    corr_mat = 1.0 - distance.squareform(corr_mat)
    sort_cols = [cols[j] for j in order]
    corr_mat = corr_mat[order][:,order[::-1]]
     
    n, m = corr_mat.shape
    
    im = ax0.matshow(corr_mat, cmap=cmap, extent=[0,m,0,n], aspect='auto')
    ax0.tick_params(top=False, bottom=True, left=True, right=False,
                    labeltop=False, labelbottom=True, labelleft=True, labelright=False)
    ax0.set_xticks(np.arange(0.5,m+0.5))
    ax0.set_xticklabels(sort_cols, rotation=90)
    ax0.set_yticks(np.arange(0.5,n+0.5))
    ax0.set_yticklabels(sort_cols)

    plt.sca(ax1)
    hierarchy.dendrogram(linkage, orientation='right', link_color_func=lambda k: '#000000')

    plt.sca(ax2)
    hierarchy.dendrogram(linkage, orientation='top', link_color_func=lambda k: '#000000')
    
    fig.colorbar(im, cax=ax3, label='Correlation')
                         
    for ax in (ax1, ax2):
        ax.spines['bottom'].set_visible(False) 
        ax.spines['top'].set_visible(False) 
        ax.spines['left'].set_visible(False) 
        ax.spines['right'].set_visible(False) 
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.set_yticklabels([])

    _watermark(fig)
    
    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()




def pvalFC_hists(plotdict, pdf, fontsize=8, colors=['#D00000','#0080FF'], nbins=50):
    
    ###  Splt this between comparison sets
    
    keys = sorted(plotdict.keys())
    n_rows = len(keys)
    
    fig, axs = plt.subplots(n_rows, 2, figsize = (8, (1.0+n_rows)*1.2), sharex='col', squeeze=False)
    
    c0, c1 = colors
    x0_min = -1.0/nbins
    x0_max = 1.0 - x0_min
     
    for row, key in enumerate(keys):
         title = key.replace(':::', ' vs ')
         df = plotdict[key]
         ax0 = axs[row, 0]
         ax1 = axs[row, 1]
         
         pvals = df['pvals'].to_list()
         FCs = df['grp1grp2_zFC'].to_list() #changed from 'grp1grp2_FC' to 'grp1grp2_zFC'
         
         if row == 0:
             ax0.set_title(f'P-value distrbutions', color=c0)    
             ax1.set_title(f'Fold change distributions', color=c1)
 
         txt0 = ax0.text(0.05, 0.8, title, transform=ax0.transAxes, fontsize=12)
         txt1 = ax1.text(0.05, 0.8, title, transform=ax1.transAxes, fontsize=12)
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
         
    
    plt.subplots_adjust(wspace=0.15, hspace=0.1, top=0.8, bottom=0.2, left=0.1, right=0.95)
    
    _watermark(fig)

    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()


def volcano(pdf, pair_name, df, q95, FClim, pthresh, colors=['#0080FF','#A0A000','#D00000'],
            split_x=True, hq_only=False, hit_labels=True, markers=None,
            lw=0.25, ls='--', lcolor='#808080', max_size=160.0, minsize=8.0,
            label_size=5.0):
    
    colors=['#0080FF','#A0A000','#D00000']
    cmap = LinearSegmentedColormap.from_list('volc', colors)
    
    group1, group2 = pair_name.split(':::')
    q951, q952 = q95
    
    if not markers:
        markers = []
   
    pthresh = -np.log2(0.01 * pthresh)
    
    high_qual = (np.array(df['zstd_grp1_q95']) != q951) & (np.array(df['zstd_grp2_q95']) != q952)
    
    has_zeros = (np.array(df['zstd_grp1_q95']) == 0.0) | (np.array(df['zstd_grp2_q95']) == 0.0)
    
    names = np.array(df.index)
    
    pvalues = np.array(df['Tpvals'])
    pvalues_q95 = np.array(df['Tpvals_q95'])
    
    nan_mask = np.isnan(pvalues)
    pvalues[nan_mask] = pvalues_q95[nan_mask]
    
    fcvalues = np.array(df['grp1grp2_zFC']) #changed from 'grp1grp2_FC' to 'grp1grp2_zFC'
    
    if hq_only:
        pvalues = pvalues[high_qual]
        fcvalues = fcvalues[high_qual]
        names = names[high_qual]
        high_qual = high_qual[high_qual]
        has_zeros = high_qual[high_qual]
     
    pos_overFC = []
    neg_overFC = []
    
    color_neg, color_low, color_pos = colors
    
    if split_x:
         xlims, neg_split, pos_split = util.get_splitxlist(df, abs(split_x))
    
    else:
         xlims = util.get_xlist(df)
         neg_split, pos_split = None, None
    
    y_label = 'p-value'
    p_label = f'{2.0**-pthresh:.2g}'
    p_label_kw = dict(xytext=(2,2), xycoords='data', textcoords='offset points', color='#808080', fontsize=5, ha='left', va='bottom')
    
    figsize=(10,8)
    xmin, xmax = xlims
    hline_kw = dict(ls=ls, alpha=0.5, lw=lw, color=lcolor)
    
    d = .8
    diag_kwargs = dict(marker=[(-1, -d), (1, d)], markersize=8,
                      linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    
    if neg_split and pos_split:
        n1, n2 = neg_split
        p1, p2 = pos_split
        width_ratios = {'width_ratios':[n1-xmin, p1-n2, xmax-p2]}
        fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=figsize, sharey=True, gridspec_kw=width_ratios)
        ax0.set_ylabel(y_label)
        ax0.set_xlim(xmin, n1)
        ax1.set_xlim(n2, p1)
        ax2.set_xlim(p2, xmax)
        ax0.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        #ax2.yaxis.tick_right()
        #ax2.yaxis.set_label_position("right")
        ax0.axhline(pthresh, **hline_kw)
        ax1.axhline(pthresh, **hline_kw)
        ax2.axhline(pthresh, **hline_kw)
        ax0.annotate(p_label, (xmin, pthresh), **p_label_kw)
        
        ax0.plot([1,1], [1,0], transform=ax0.transAxes, **diag_kwargs)
        ax1.plot([0,0], [1,0], [1,1], [1,0], transform=ax1.transAxes, **diag_kwargs)
        ax2.plot([0,0], [0,1], transform=ax2.transAxes, **diag_kwargs)
        
        if xmin < -FClim <= n1:
            ax0.axvline(-FClim,   **hline_kw)
        elif n2 <= -FClim < 0.0:
            ax1.axvline(-FClim,   **hline_kw)
        
        if 0.0 < FClim <= p1:
            ax1.axvline(FClim,   **hline_kw)
        elif p2 <= FClim < xmax:
            ax2.axvline(FClim,   **hline_kw)
        
    elif neg_split:
        n1, n2 = neg_split
        width_ratios = {'width_ratios':[n1-xmin, xmax-n2]}
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=figsize, sharey=True, gridspec_kw=width_ratios)
        ax0.set_ylabel(y_label)
        ax0.set_xlim(xmin, n1)
        ax1.set_xlim(n2, xmax)
        ax2 = None
        ax0.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax0.axhline(pthresh, **hline_kw)
        ax1.axhline(pthresh, **hline_kw)
        ax0.annotate(p_label, (xmin, pthresh), **p_label_kw)

        ax0.plot([1,1], [1,0], transform=ax0.transAxes, **diag_kwargs)
        ax1.plot([0,0], [1,0], transform=ax1.transAxes, **diag_kwargs)
        
        if xmin < -FClim <= n1:
            ax0.axvline(-FClim,   **hline_kw)
        elif n2 <= -FClim < xmax:
            ax1.axvline(-FClim,   **hline_kw)
        
    elif pos_split:
        p1, p2 = pos_split
        width_ratios = {'width_ratios':[p1-xmin, xmax-p2]}
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, sharey=True, gridspec_kw=width_ratios)
        ax0 = None
        ax1.set_ylabel(y_label)
        ax1.set_xlim(xmin, p1)
        ax2.set_xlim(p2, xmax)
        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        #ax2.yaxis.tick_right()
        #ax2.yaxis.set_label_position("right")
        ax1.axhline(pthresh, **hline_kw)
        ax2.axhline(pthresh, **hline_kw)
        ax1.annotate(p_label, (xmin, pthresh), **p_label_kw)

        ax1.plot([1,1], [1,0], transform=ax1.transAxes, **diag_kwargs)
        ax2.plot([0,0], [0,1], transform=ax2.transAxes, **diag_kwargs)
        
        if xmin < FClim <= p1:
            ax1.axvline(FClim,   **hline_kw)
        elif p2 < FClim <= xmax:
            ax2.axvline(FClim,   **hline_kw)
        
    else:
        fig, ax1 = plt.subplots(figsize=figsize)
        ax0, ax2 = None, None
        ax1.set_ylabel(y_label)
        ax1.set_xlim(xmin, xmax)
        ax1.axhline(pthresh, **hline_kw)
        ax1.axvline(FClim,   **hline_kw)
        ax1.axvline(-FClim,  **hline_kw)
        ax1.annotate(p_label, (xmin, pthresh), **p_label_kw)
        
    fig.subplots_adjust(wspace=0.05)
    fig.set_clip_on(False)
    
    ax1.set_title(pair_name.replace(':::',' vs '), fontsize=18)
    ax1.set_xlabel('log2 fold change')
    
    pnorm = np.array(pvalues)
    pnorm /= pnorm.max() # 0..1
    
    fcnorm = np.abs(fcvalues)
    fcnorm /= fcnorm.max() # 0..1
    
    weights = pnorm * fcnorm # -1..1
    weights /= np.abs(weights).max()
    weights = (1.0 + weights) / 2.0 # 0..1
    
    sizes = max_size * pnorm * pnorm * fcnorm * fcnorm # 0 .. max_size
    sizes = np.clip(sizes, minsize, max_size)
    
    
    scatter_kw = dict(alpha=0.4, edgecolors='none', clip_on=False)
    hit_kw = dict(alpha=1.0, clip_on=False, edgecolor='k', linewidth=0.5, zorder=10)
   
    # middle bounds
    m0 = xmin
    m1 = xmax
 
    if ax0:
        selection = fcvalues < neg_split[0]
        
        if np.any(selection):
            colors = [cmap(w) for w in weights[selection]]
            ax0.scatter(fcvalues[selection], pvalues[selection], s=sizes[selection], c=colors,**scatter_kw)
        
        m0 = neg_split[1]
        
    if ax2:
        selection = fcvalues > pos_split[1]
        if np.any(selection):
            colors = [cmap(w) for w in weights[selection]]
            ax2.scatter(fcvalues[selection], pvalues[selection], s=sizes[selection], c=colors, **scatter_kw)
        
        m1 = pos_split[0]
    
    selection = (fcvalues >= m0) & (fcvalues <= m1) & high_qual
    
    if np.any(selection):
        colors = [cmap(w) for w in weights[selection]]
        ax1.scatter(fcvalues[selection], pvalues[selection], s=sizes[selection], c=colors, **scatter_kw)

    marker_text = []

    for i, name in enumerate(names):
        size = sizes[i]
        low_qual = not high_qual[i]
        color =  color=cmap(weights[i])
        ds = 0.5 * np.sqrt(size)
        x = fcvalues[i]
        y = pvalues[i]
        
        if hq_only and low_qual:
            continue
        
        if ax0 and (x <= neg_split[0]):
            ax = ax0
        elif ax2 and (x >= pos_split[1]):
            ax = ax2
        else:
            ax = ax1  
 
        if name in markers:
            marker_text.append((x, y, ds, name))
            ax.scatter(x, y, s=size, color='k', edgecolor='k', clip_on=False, zorder=11)
            continue
        
        if has_zeros[i]:
          name = name + '*'
        
        if x >= FClim and y >= pthresh:
            if low_qual:
                 ax.scatter(x, y, s=size, marker='h', color=color, **hit_kw)
            
            else:
                 pos_overFC.append(name)
                 ax.scatter(x, y, s=size, color=color, **hit_kw)
            
            if hit_labels:
                for a in (ax0, ax1, ax2):
                  if not a:
                    continue
 
                  txt = a.annotate(name, xy=(x, y), xytext=(ds, ds), fontsize=label_size, textcoords='offset points', clip_on=False, zorder=11)
                  txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF80')])
             
            continue
             
        if x <= -FClim and y >= pthresh:
            if low_qual:
                 ax.scatter(x, y, s=size, marker='h', color=color, **hit_kw)
            else:     
                 neg_overFC.append(name)
                 ax.scatter(x, y, s=size, color=color, **hit_kw)
                 
            if hit_labels:
                for a in (ax0, ax1, ax2):
                  if not a:
                    continue
                  txt = a.annotate(name, xy=(x, y), xytext=(-ds, ds), fontsize=label_size, textcoords='offset points', ha='right', clip_on=False, zorder=11)
                  txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF80')])
            
            continue
 
        if low_qual:
            ax.scatter(x, y, s=size, marker='v', color=color, zorder=12) # color=color_low
            continue
        
    for x, y, ds, name in marker_text:
         txt = ax.annotate(name, xy=(x, y), xytext=(ds, ds), textcoords='offset points',
                          fontsize=label_size, clip_on=False)
         txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF80')])
    
    ax = ax0 if ax0 else ax1
    ax.text(0.05, 0.05, f'Greater in {group2}', color=color_neg, ha='left', va='top', transform=ax.transAxes, fontsize=10, clip_on=False)
    
    ax.set_ylim(0.0)
    pmax = math.ceil(np.nanmax(pvalues))
    p = int(pmax) 
    y_ticks = np.logspace(-p, 0, p+1, base=10.0)
    y_ticks = y_ticks[y_ticks > 2.0 ** (-pmax)]
    
    ax.yaxis.tick_left()
    ax.set_yticks(-np.log2(y_ticks))
    ax.set_yticklabels(y_ticks)
    
    ax = ax2 if ax2 else ax1
    ax.text(0.95, 0.05, f'Greater in {group1}', color=color_pos, ha='right', va='top', transform=ax.transAxes, fontsize=10, clip_on=False)
    
    ax1.set_facecolor('none') # Stops covering part of text
    
    if ax0:
        ax1.tick_params('y', length=0, color='w')
        ax0.set_facecolor('none')
    if ax2:
        ax2.tick_params('y', length=0, color='w')
        ax2.set_facecolor('none')

    _watermark(fig)

    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()
    

