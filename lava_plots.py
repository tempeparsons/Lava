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
            mu = np.mean(data)
            sig = np.std(data)
            t = f'n={n:,}\n$\mu$={mu:.2}\n$\sigma$={sig:.2}'
            
            txt = ax.text(0.95, 0.9, t, ha='right', va='top', color=color, transform=ax.transAxes, fontsize=8)
            txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
            
            if (i // nx) == ny-1:
                ax.set_xlabel('Data value')
            
            if (i % nx) == 0:
                ax.set_ylabel('Bin count')
            
        else:
            ax.set_visible(False)

    ax.set_ylim(0.0)
    """
    df.hist(column=value_cols, 
            alpha=0.5, 
            bins=50, 
            grid=True, 
            sharex=True, sharey=True, 
            xlabelsize=6, ylabelsize=6,
            figsize=(9,9))
    """
    plt.subplots_adjust(hspace=0.05, wspace=0.05, top=0.95, bottom=0.1, left=0.1, right=0.95) 
    
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
 
                #txt = ax.text(0.05, 0.9, f'{cols[i]} vs {cols[j]}', transform=ax.transAxes, fontsize=12)
                #txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
 
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

    """
    nrows = df.shape[1] // ncols + (df.shape[1] % ncols > 0)
    cols = value_cols
    fig, axs = plt.subplots(nrows, ncols, figsize=(9,9), sharex=True, sharey=True)
    for i in range(len(cols)-1):
        axs.flat[i].scatter(df[cols[i]], df[cols[i+1]], s=5, alpha=0.6)
        txt = axs.flat[i].set_title(f'{cols[i]} vs {cols[i+1]}', fontsize=8)
        txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
        ymin, ymax = axs.flat[i].get_ylim()
        axs.flat[i].axline([ymax*0.5,ymax*0.5], slope=1, c='r', lw=0.5)    
    # switch off unused axes
    for i in range(len(value_cols)-1, nrows*ncols):
        axs.flat[i].set_visible(False)
    """
        
    plt.subplots_adjust(wspace=0.0, hspace=0.0, top=0.95, bottom=0.05, left=0.05, right=0.95)
    
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
    
    if pdf:
      pdf.savefig(dpi=DPI)
    else:
      plt.show()
    
    plt.close()




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
         
    
    plt.subplots_adjust(wspace=0.15, hspace=0.1, top=0.95, bottom=0.05, left=0.1, right=0.95)

    if pdf:
      pdf.savefig(dpi=DPI)
    else:
      plt.show()
    
    plt.close()


def volcano(pdf, pair_name, df, q95, FClim, pthresh, colors=['#0080FF','#A0A000','#D00000'],
            hq_only=False, hit_labels=True, print_hits=True, markers=None,
            lw=0.25, ls='--', lcolor='#808080', max_size=160.0, minsize=8.0):
    
    colors=['#0080FF','#A0A000','#D00000']
    cmap = LinearSegmentedColormap.from_list('volc', colors)
    
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
    
    color_neg, color_low, color_pos = colors
    xlims = util.get_xlist(df)
 
    fig, ax = plt.subplots(figsize = (10,8), sharey=True)
 
    ax.set_xlim(xlims[0], xlims[1])

    ax.axvline(FClim,   ls=ls, alpha=0.5, lw=lw, color=lcolor)
    ax.axvline(-FClim,  ls=ls, alpha=0.5, lw=lw, color=lcolor)
    ax.axhline(pthresh, ls=ls, alpha=0.5, lw=lw, color=lcolor)

    ax.set_title(pair_name, fontsize=18)
    #ax.set_ylabel('-log2 transformed p-value')
    ax.set_ylabel('p-value')
    ax.set_xlabel('log2 fold change')
    
    pnorm = plotlist[1]
    pnorm /= pnorm.max() # 0..1
    
    fcnorm = np.abs(plotlist[2])
    fcnorm /= fcnorm.max() # 0..1
    
    weights = pnorm * fcnorm # -1..1
    sizes = max_size * pnorm * pnorm * fcnorm * fcnorm # 0 .. max_size
    sizes = np.clip(sizes, minsize, max_size)
    
    ax.scatter(plotlist[1], plotlist[2], s=sizes, c=weights, cmap=cmap, norm=Normalize(-1.0, 1.0), alpha=0.4, edgecolors='none')

    marker_text = []

    for i, name in enumerate(plotlist[0]):
        size = sizes[i]
        ds = 0.5 * np.sqrt(size)
        
        if name in markers:
            marker_text.append((plotlist[1][i], plotlist[2][i], name))
            ax.scatter(plotlist[1][i], plotlist[2][i], s=3, color='k', edgecolor='k')

        if plotlist[1][i] >= FClim and plotlist[2][i] >= pthresh:
            pos_overFC.append(name)
            ax.scatter(plotlist[1][i], plotlist[2][i], s=size, color=cmap(weights[i]), edgecolor='k', linewidth=0.5)
            if hit_labels:
                txt = ax.annotate(name, xy=(plotlist[1][i], plotlist[2][i]), xytext=(ds, ds), fontsize=6, textcoords='offset points')
                txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])

        if plotlist[1][i] <= -FClim and plotlist[2][i] >= pthresh:
            neg_overFC.append(name)
            ax.scatter(plotlist[1][i], plotlist[2][i], s=size, color=cmap(weights[i]), edgecolor='k', linewidth=0.5)
            if hit_labels:
                txt = ax.annotate(name, xy=(plotlist[1][i], plotlist[2][i]), xytext=(ds, ds), fontsize=6, textcoords='offset points')
                txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])

        if not hq_only:
            if plotlist[3][i] == q951 or plotlist[4][i] == q952:
                ax.scatter(plotlist[1][i], plotlist[2][i], s=6, color=color_low)

    for tx in marker_text:
            ax.annotate(tx[2], xy=(tx[0], tx[1]),  xytext=(tx[0], tx[1]),
                    arrowprops=dict(arrowstyle='-', fc="k", ec="k", lw=0.5, relpos=(0.25, 0.5)),
                    bbox=dict(pad=-2, facecolor="none", edgecolor="none"),
                    ha="left", va="center", size=6)
 
    #legend_elements = [Line2D([0], [0], marker='o', color='w', label=f'Greater in {group1}',markerfacecolor=color_pos, markersize=8, linestyle=''),
    #                   Line2D([0], [0], marker='o', color='w', label=f'Greater in {group2}',markerfacecolor=color_neg, markersize=8, linestyle='')]
    #ax.legend(handles=legend_elements)
    
    ax.text(0.95, 0.05, f'Greater in {group1}', color=color_pos, ha='right', va='top', transform=ax.transAxes, fontsize=10)
    ax.text(0.05, 0.05, f'Greater in {group2}', color=color_neg, ha='left', va='top', transform=ax.transAxes, fontsize=10)
    
    pmax = math.ceil(np.max(plotlist[2]))
    
    ax.set_ylim(0.0)
    y_ticks = np.logspace(-9, 0, 10, base=10.0)
    y_ticks = y_ticks[y_ticks > 2.0 ** (-pmax)]
    
    ax.set_yticks(-np.log2(y_ticks))
    ax.set_yticklabels(y_ticks)
    
    if print_hits:
        print('greater in grp1:', pos_overFC)
        print('greater in grp2:', neg_overFC)

    if pdf:
      pdf.savefig(dpi=DPI)
    else:
      plt.show()
    
    plt.close()
    

