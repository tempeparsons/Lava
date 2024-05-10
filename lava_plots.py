import math
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.patheffects as PathEffects
from matplotlib.lines import Line2D
from matplotlib.text import Text
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


def _nice_tick_labels(labels, texp=1e-4):
    
    new = []
    for t in labels:
        if t > texp:
            new.append(f'{t}')
        else:
            pwr = np.log10(t) 
            new.append(f'$10^{{{pwr}}}$')
    
    
    return new
    
    
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
    
    maxw = max([len(x) for x in list(value_cols)])
    fontsize = min(8, 80.0/maxw)
    
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,8))
    df.boxplot(column=value_cols, ax=ax1, grid=False, rot=45, fontsize=fontsize, boxprops=box_style,
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
        ax2.scatter(x_vals, y_vals, color=color, alpha=0.3, s=4)
 
    x_ticks = np.arange(0, nc)   
    ax2.set_xlim(-0.5, nc-0.5)  
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_ticks+1)
    ax2.set_xlabel('Sample column', fontsize=9)
    ax2.set_ylabel('$log_2$(abundance)', fontsize=9)
    ax1.set_ylabel('$log_2$(abundance)', fontsize=9)
    ax2.tick_params('x', top=True)
    
    plt.subplots_adjust(hspace=0.25, top=0.90, bottom=0.1, left=0.10, right=0.95)    
    _watermark(fig)
    
    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()
    
    
def plot_dual_pca(df, value_cols, norm_cols, group_cols, pdf, colors=['#D00000','#A0A000','#0080FF'], size=10):
    
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    
    groups = sorted(group_cols)
    n_groups = len(groups)
    col_groups = {}
    for group, cols in group_cols.items():
        for col in cols:
            col_groups[col] = group
    
    data = np.array(df[value_cols])
    zdata = np.array(df[norm_cols])
    
    n_prots, n_samples = data.shape
    
    # Remove all blank
    s_means = np.nanmean(data, axis=0)
    data = np.nan_to_num(data)
    data = data[data.sum(axis=1) > 0.0]
    sz_means = np.nanmean(zdata, axis=0)
    zdata = np.nan_to_num(zdata)
    zdata = zdata[zdata.sum(axis=1) > 0.0]
    
    for i in range(n_samples):
        z = data[:,i] == 0.0
        data[z, i] = s_means[i]
        z = zdata[:,i] == 0.0
        zdata[z, i] = sz_means[i]
    
    scaler = StandardScaler()
    data = data.T  # expt_samples, proteins
    data = scaler.fit_transform(data)

    zdata = zdata.T  # expt_samples, proteins
    zdata = scaler.fit_transform(zdata)

    cmap = LinearSegmentedColormap.from_list('group_colors', colors)
    
    grp_colors = [cmap(x) for x in np.linspace(0.0, 1.0, n_groups)]
    grp_colors = dict(zip(groups, grp_colors))
    
    sample_colors = [grp_colors[col_groups[col]] for col in value_cols]
    
    pca = PCA(n_components=2)
    data_trans = pca.fit_transform(data)
    e_var_x, e_var_y = pca.explained_variance_ratio_
    
    pca2 = PCA(n_components=2)
    zdata_trans = pca2.fit_transform(zdata)
    e_var_zx, e_var_zy = pca2.explained_variance_ratio_

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    
    x_vals, y_vals = data_trans.T
    x_zvals, y_zvals = zdata_trans.T
    
    ax1.set_title('Input values PCA')
    ax1.scatter(x_vals, y_vals, s=size, color=sample_colors)
    ax1.set_xlabel(f'PC1 ({e_var_x*100.0:.2f}%)')
    ax1.set_ylabel(f'PC2 ({e_var_y*100.0:.2f}%)')

    ax2.set_title('Background adjusted PCA')
    ax2.scatter(x_zvals, y_zvals, s=size, color=sample_colors)
    ax2.set_xlabel(f'PC1 ({e_var_zx*100.0:.2f}%)')
    ax2.set_ylabel(f'PC2 ({e_var_zy*100.0:.2f}%)')
    
    for group in groups:
        ax1.scatter([], [], s=size, color=grp_colors[group], label=group)
    
    for i, col in enumerate(value_cols):
        txt = ax1.annotate(col, xy=(x_vals[i], y_vals[i]), xytext=(5, 5), fontsize=5, alpha=0.5,
                           ha='left', va='top', textcoords='offset points', clip_on=False)
        txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
        txt = ax2.annotate(col, xy=(x_zvals[i], y_zvals[i]), xytext=(5, 5), fontsize=5, alpha=0.5,
                           ha='left', va='top', textcoords='offset points', clip_on=False)
        txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])

    plt.subplots_adjust(wspace=0.25, top=0.90, bottom=0.15, left=0.10, right=0.95)    
    
    ax1.legend(fontsize=5) # , bbox_to_anchor=(1.04, 1), loc='upper left', )

    _watermark(fig)
    
    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()
    

    
def plot_pca(df, value_cols, group_cols, pdf, colors=['#D00000','#A0A000','#0080FF'], size=10):
    
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    
    groups = sorted(group_cols)
    n_groups = len(groups)
    col_groups = {}
    for group, cols in group_cols.items():
        for col in cols:
            col_groups[col] = group
    
    data = np.array(df[value_cols])
    n_prots, n_samples = data.shape
    
    # Remove all blank
    s_means = np.nanmean(data, axis=0)
    data = np.nan_to_num(data)
    data = data[data.sum(axis=1) > 0.0]
    #data = data[np.count_nonzero(data, axis=1) == n_samples]
    
    for i in range(n_samples):
        z = data[:,i] == 0.0
        data[z, i] = s_means[i]
    
    scaler = StandardScaler()
    data = data.T  # expt_samples, proteins
    data = scaler.fit_transform(data)

    cmap = LinearSegmentedColormap.from_list('group_colors', colors)
    
    grp_colors = [cmap(x) for x in np.linspace(0.0, 1.0, n_groups)]
    grp_colors = dict(zip(groups, grp_colors))
    
    cols = [x[6:] if x.startswith('znorm_') else x for x in value_cols]
    sample_colors = [grp_colors[col_groups[col]] for col in cols]
    
    pca = PCA(n_components=2)

    data_trans = pca.fit_transform(data)
    
    e_var_x, e_var_y = pca.explained_variance_ratio_

    fig, ax1 = plt.subplots(1, 1, figsize=(7, 5))
    
    x_vals, y_vals = data_trans.T
    
    ax1.set_title('Sample PCA')
    ax1.scatter(x_vals, y_vals, s=size, color=sample_colors)
    ax1.set_xlabel(f'PC1 ({e_var_x*100.0:.2f}%)')
    ax1.set_ylabel(f'PC2 ({e_var_y*100.0:.2f}%)')
    
    for group in groups:
        ax1.scatter([], [], s=size, color=grp_colors[group], label=group)
    
    for i, col in enumerate(cols):
        txt = ax1.annotate(col, xy=(x_vals[i], y_vals[i]), xytext=(5, 5), fontsize=7, alpha=0.5,
                           ha='left', va='top', textcoords='offset points', clip_on=False)
        txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])

    plt.subplots_adjust(hspace=0.25, top=0.90, bottom=0.1, left=0.10, right=0.8)    
    
    ax1.legend(fontsize=7, bbox_to_anchor=(1.04, 1), loc='upper left', )

    _watermark(fig)
    
    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()
    
    
    
def plot_norm(df, value_cols, norm_cols, pdf, colors=['#D00000','#A0A000','#0080FF'], bin_size=0.5, alpha=0.7):

    cmap = LinearSegmentedColormap.from_list('main', colors)

    fig, (ax1, ax2, ax3) = plt.subplots(1,3, sharey=True, figsize=(8, 4), gridspec_kw={'width_ratios': [3, 3, 1]})
    
    n = len(value_cols)

    all_vals = np.array(df[value_cols])
    min_val = int(10.0 * np.nanmin(all_vals)) * 0.1
    max_val = int(10.0 * np.nanmax(all_vals)) * 0.1
    bins = np.arange(min_val, max_val, bin_size)

    all_zvals = np.array(df[norm_cols])
    min_zval = int(10.0 * np.nanmin(all_zvals)) * 0.1
    max_zval = int(10.0 * np.nanmax(all_zvals)) * 0.1
    zbins = np.arange(min_zval, max_zval, bin_size)
    
    for i in range(n):
       color = cmap(i/(n-1.0))
    
       data = df[value_cols[i]]
       hist, edges = np.histogram(data[~np.isnan(data)], bins=bins, density=True)
       x_vals =  0.5 * (edges[:-1] + edges[1:])
       ax1.plot(x_vals, hist, color=color, linewidth=0.50, alpha=alpha)

       data = df[norm_cols[i]]
       hist, edges = np.histogram(data[~np.isnan(data)], bins=zbins, density=True)
       x_vals =  0.5 * (edges[:-1] + edges[1:])
       ax2.plot(x_vals, hist, color=color, linewidth=0.5, alpha=alpha)
       
       ax3.plot([], color=color, linewidth=0.5, alpha=alpha, label=value_cols[i])
      
    ax1.set_title('Input values')
    ax1.set_xlabel('$log_2$(abundance)', fontsize=9)
    ax1.set_ylabel('Probability density', fontsize=9)
    
    ax2.set_title('Normalized values')
    ax2.set_xlabel('Probability density', fontsize=9)
    ax2.set_xlabel('Norm. $log_2$(abundance)', fontsize=9)
    
    ax3.set_axis_off()
    ax3.legend(fontsize=5)
    
    plt.subplots_adjust(wspace=0.05, top=0.92, bottom=0.12, left=0.10, right=0.95)    
    _watermark(fig)
    
    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()
    
    
def histograms(df, value_cols, pdf, colors=['#D00000','#A0A000','#0080FF'], bin_size=0.5):

    cols = list(df[value_cols])
    n_cols = len(cols)
    
    nx = int(math.ceil(math.sqrt(n_cols)))
    ny = int(math.ceil(n_cols/nx))
    
    fig, axarr = plt.subplots(ny, nx, squeeze=False, figsize=(9, 7), sharex=True, sharey=True)
    
    all_vals = np.array(df[value_cols])
    min_val = int(10.0 * np.nanmin(all_vals)) * 0.1
    max_val = int(10.0 * np.nanmax(all_vals) + 1.0) * 0.1
    bins = np.arange(min_val, max_val, bin_size)
     
    for i in range(nx*ny):
        ax = axarr.flat[i]
        ax.tick_params(axis='both', which='major', labelsize=8)
        
        if i < n_cols:
            color = colors[i % len(colors)]
            data = df[value_cols[i]]
            data = data[~np.isnan(data)]
            hist, edges = np.histogram(data[~np.isnan(data)], bins=bins)
            
            x_vals =  0.5 * (edges[:-1] + edges[1:])
            ax.plot(x_vals, hist, color=color, linewidth=1.0)
            ax.fill_between(x_vals, hist, 0, color=color, alpha=0.5)
            
            fontsize = min(12, 120.0/len(cols[i]))
            
            txt = ax.text(0.05, 0.9, cols[i], transform=ax.transAxes, fontsize=fontsize)
            txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
            
            n = len(data)
            #mu = np.mean(data)
            #sig = np.std(data)
            t = f'n={n:,}' #'\n$\mu$={mu:.2}\n$\sigma$={sig:.2}'
            
            txt = ax.text(0.95, 0.9, t, ha='right', va='top', color=color, transform=ax.transAxes, fontsize=8)
            txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')])
            
            if (i // nx) == ny-1:
                ax.set_xlabel('$log_2$(abundance)', fontsize=9)
            
            if (i % nx) == 0 and ((ny < 3) or ((i // ny) == (ny // 2))):
                ax.set_ylabel(f'Bin count (width {bin_size})', fontsize=9)
            
            ax.grid(visible=True, which='major', axis='x', color='#808080', linewidth=0.5, linestyle='--', alpha=0.5, zorder=-1)
            
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
        
    
def xy_plots(df, value_cols, ncols, pdf, colors=['#00FFFF','#0000FF','#000000','#A00000','#AAAA00'], figsize=10.0, bg_color='#B0B0B0'):
    
    main_cmap = LinearSegmentedColormap.from_list('main', colors) # bad to good
    
    cols = list(df[value_cols])
    n = len(cols)    
    
    fig, axarr = plt.subplots(n, n, squeeze=False, figsize=(figsize,figsize), sharex=True, sharey=True)
    
    gridsize = int(400//n)
    fontsize = min(1e3/(figsize * n), 14)
    
    max_label = max([len(x) for x in list(value_cols)])
    label_fontsize = min(fontsize, (10.0*fontsize)/max_label)
    
    data = np.array([df[col] for col in cols])
    data[~(data > 0)] =  np.nanmean(data) # Zeros and nans
    corr_mat = distance.pdist(data, metric='correlation') # Correlation 'distance" : 1.0 - rho
    linkage = hierarchy.ward(corr_mat)
    order = hierarchy.leaves_list(linkage)
    
    corr_mat = 1.0 - distance.squareform(distance.pdist(data, metric='correlation'))
    min_corr = corr_mat.min()
    max_corr = corr_mat.max()
     
    for row, i in enumerate(order):
        data1 = np.array(df[cols[i]])
        
        for col, j in enumerate(order):
            ax = axarr[n-row-1,col]

            data2 = np.array(df[cols[j]])
            valid = (data1 > 0) & (data2 > 0) # no nans or zeros 
            
            x_vals = data1[valid]
            y_vals = data2[valid]
            
            rho = corr_mat[i,j]
            goodness = (rho + 1.0) / 2.0 # (rho-min_corr)/(max_corr-min_corr)
            
            cmap = LinearSegmentedColormap.from_list('{i}+{j}', [main_cmap(goodness), '#000000'])
            ax.hexbin(x_vals, y_vals, cmap=cmap, gridsize=gridsize, mincnt=1, edgecolors='none')

            if row > col:
                txt = ax.text(0.5, 0.5, f'{rho:.2}', ha='center', va='center', transform=ax.transAxes, fontsize=1.5*fontsize)
                txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF')]) 
            
            ax.set_xticks([])
            ax.set_yticks([])
            
            if row == 0:
                ax.set_xlabel(f'{cols[j]}', fontsize=label_fontsize)
                
            elif row == (n-1):    
                ax.xaxis.set_label_position('top') 
                ax.set_xlabel(f'{cols[j]}', fontsize=label_fontsize)
  
            if col == 0:
                ax.set_ylabel(f'{cols[i]}', fontsize=label_fontsize)
                
            elif col == (n-1):
                ax.yaxis.set_label_position('right')     
                ax.set_ylabel(f'{cols[i]}', fontsize=label_fontsize)

    plt.subplots_adjust(wspace=0.0, hspace=0.0, top=0.95, bottom=0.05, left=0.05, right=0.95)
    _watermark(fig)
    
    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()


def correlation_plot(df, value_cols, pdf, colors=['#00FFFF','#0000FF','#000000','#A00000','#FFFF00'],
                     figsize=9.0, bg_color='#B0B0B0', title='Sample correlation matrix'):
    
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
    ax2.set_title(title)
    # Colorbar
    ax3 = fig.add_axes([0.87, 0.85, 0.05, 0.10])
    
    data = np.array([df[col] for col in cols])
    data[~(data > 0)] =  np.nanmean(data) # Zeros and nans
    
    corr_mat = distance.pdist(data, metric='correlation') # Correlation 'distance" : 1.0 - rho
    
    linkage = hierarchy.ward(corr_mat)
    order = hierarchy.leaves_list(linkage)
    
    corr_mat = 1.0 - distance.squareform(corr_mat)
    sort_cols = [cols[j] for j in order]
    corr_mat = corr_mat[order[::-1]][:,order]
     
    n, m = corr_mat.shape
   
    max_label = max([len(x) for x in list(value_cols)])
    label_fontsize = min(12, 120.0/max_label)
    
    im = ax0.matshow(corr_mat, cmap=cmap, extent=[0,m,0,n], aspect='auto', vmin=-1.0, vmax=1.0)
    ax0.tick_params(top=False, bottom=True, left=True, right=False,
                    labeltop=False, labelbottom=True, labelleft=True, labelright=False)
    ax0.set_xticks(np.arange(0.5,m+0.5))
    ax0.set_xticklabels(sort_cols, rotation=90, fontsize=label_fontsize)
    ax0.set_yticks(np.arange(0.5,n+0.5))
    ax0.set_yticklabels(sort_cols, fontsize=label_fontsize)

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
         
         pvals = np.array(df['pvals'])
         FCs = np.array(df['grp1grp2_FC'])
         
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
         
    
    plt.subplots_adjust(wspace=0.15, hspace=0.15, top=0.8, bottom=0.2, left=0.1, right=0.95)
    
    _watermark(fig)

    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()


def pp_plot(pdf, pair1, df1, pair2, df2, p_thresh, fc_thresh, min_peps, colors=['#0080FF','#A0A000','#D00000'],
            max_size=160.0, min_size=8.0, label_size=5.0):

    # ? Add option to have these plots
    # Table output: fcs and pvals, max(pvalue) > thresh, max(p_value) > thresh (different colours), n_obs 1,2,3,4, n_pep, samce accession, labels, extra cols as volc.
    
    a1, a2 = pair1
    b1, b2 = pair2
    nlp = - np.log2(p_thresh)
    
    cmap = LinearSegmentedColormap.from_list('pp', colors)
    grey = '#80808080'
    lcolor = '#808080'
    
    hline_kw = dict(ls='--', alpha=0.5, lw=0.5, color=lcolor)
    index1 = df1.index
    index2 = df2.index
    
    common = set(index1) & set(index2)
    not1 = set(index2) - set(index1)
    not2 = set(index1) - set(index2)
    
    pvalues1 = np.array(df1.loc[common, 'Tpvals'])
    fcvalues1 = np.array(df1.loc[common, 'grp1grp2_FC'])
    
    pvalues2 = np.array(df2.loc[common, 'Tpvals'])
    fcvalues2 = np.array(df2.loc[common, 'grp1grp2_FC'])
    
    names = np.array(df1.loc[common, 'labels'])
    
    has_zeros = (np.array(df1.loc[common, 'nobs_grp1']) == 0) | \
                (np.array(df1.loc[common, 'nobs_grp2']) == 0) | \
                (np.array(df2.loc[common, 'nobs_grp1']) == 0) | \
                (np.array(df2.loc[common, 'nobs_grp2']) == 0)
    
    low_qual =  (np.array(df1.loc[common, 'nobs_grp1']) < 2) | \
                (np.array(df1.loc[common, 'nobs_grp2']) < 2) | \
                (np.array(df2.loc[common, 'nobs_grp1']) < 2) | \
                (np.array(df2.loc[common, 'nobs_grp2']) < 2)
   
    if 'npeps' in df1:
        npeps = np.array(df1.loc[common, 'npeps'])
    else:
        npeps = np.full(len(names), min_peps)
    
    dual = np.ones(len(names), bool)
    
    if not1: # Not in first comparison
        n = len(not1)
        pvalues1 = np.concatenate([pvalues1, np.zeros(n)])
        fcvalues1 = np.concatenate([fcvalues1, np.zeros(n)])
        pvalues2 = np.concatenate([pvalues2, np.array(df2.loc[not1, 'Tpvals'])])
        fcvalues2 = np.concatenate([fcvalues2,np.array(df2.loc[not1, 'grp1grp2_FC'])])
        names = np.concatenate([names,np.array(df2.loc[not1, 'labels'])])
        has_zeros = np.concatenate([has_zeros, np.zeros(n, bool)])
        low_qual = np.concatenate([low_qual, np.zeros(n, bool)])
        dual = np.concatenate([dual, np.zeros(n, bool)])
        if 'npeps' in df1:
            npeps = np.concatenate([npeps, np.array(df2.loc[not1, 'npeps'])])
     
    if not2: # Not in second comparison
        n = len(not2)
        pvalues2 = np.concatenate([pvalues2, np.zeros(n)])
        fcvalues2 = np.concatenate([fcvalues2, np.zeros(n)])
        pvalues1 = np.concatenate([pvalues1, np.array(df1.loc[not2, 'Tpvals'])])
        fcvalues1 = np.concatenate([fcvalues1, np.array(df1.loc[not2, 'grp1grp2_FC'])])
        names = np.concatenate([names, np.array(df1.loc[not2, 'labels'])])
        has_zeros = np.concatenate([has_zeros, np.zeros(n, bool)])
        low_qual = np.concatenate([low_qual, np.zeros(n, bool)])
        dual = np.concatenate([dual, np.zeros(n, bool)])

        if 'npeps' in df1:
            npeps = np.concatenate([npeps, np.array(df1.loc[not2, 'npeps'])])
    
    p_best = np.maximum(pvalues1, pvalues2)
    f_worst = np.where(np.abs(fcvalues1) < np.abs(fcvalues2), fcvalues1, fcvalues2)
    
    is_hit = ((pvalues1 >= nlp) & (np.abs(fcvalues1) >= fc_thresh)) | ((pvalues2 >= nlp) & (np.abs(fcvalues2) >= fc_thresh)) 
    
    p_best /= p_best.max() # 0..1
    f_worst /= np.abs(f_worst).max() # 0..1
    
    weights = f_worst # -1..1
    weights /= np.abs(weights).max()

    sizes = weights * weights
    sizes = min_size + (max_size-min_size) * sizes
    
    weights = (1.0 + weights) / 2.0 # 0..1
    
    p_min = 0.0
    p_max = math.ceil(max(pvalues1.max(), pvalues2.max()))
    
    s0 = (npeps < min_peps) & dual 
    s1 = ~is_hit & ~s0 & dual
    s2 = is_hit & ~s0 & dual
    s3 = is_hit & ~s0 & ~dual
    
    pvalues1[fcvalues1 < 0.0] *= -1.0
    pvalues2[fcvalues2 < 0.0] *= -1.0

    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_title(f'[{a1} vs {a2}] vs [{b1} vs {b2}]', fontsize=14)
    
    ax.scatter(pvalues1[s0], pvalues2[s0], color=lcolor, alpha=0.3, clip_on=False, s=sizes[s0], zorder=1)
    ax.scatter(pvalues1[s1], pvalues2[s1], c=[cmap(w) for w in weights[s1]], alpha=0.3, clip_on=False, s=sizes[s1], zorder=2)
    ax.scatter(pvalues1[s2], pvalues2[s2], c=[cmap(w) for w in weights[s2]], alpha=1.0,  edgecolor='k', linewidth=0.5, clip_on=False, s=sizes[s2], zorder=4)
    ax.scatter(pvalues1[s3], pvalues2[s3], c=[cmap(w) for w in weights[s3]], alpha=1.0,  edgecolor='#800080', linewidth=0.5, clip_on=False, s=sizes[s3], zorder=3)
    
    hit_idx = np.nonzero(is_hit)[0]
    for i in hit_idx:
        ds = 0.5 * np.sqrt(sizes[i])
        name = names[i]
        text_alpha = 1.0
        low_pep = s0[i]
        x = pvalues1[i]
        y = pvalues2[i]
 
        if not isinstance(name, str):
            name = '?'
 
        if low_pep:
            text_color = lcolor
            text_alpha = 0.5
            zorder = 1
        
        if not dual[i]:
            name = name + '$^{!}$'
            text_color = '#800080'
            text_alpha = 0.5
            zorder = 2
          
        elif has_zeros[i]:
            name = name + '$^{0}$'
            text_color = '#400040'
            zorder = 3
 
        elif low_qual[i]:
            name = name + '$^{1}$'
            text_color = '#400040'
            zorder = 3
 
        else:
            text_color = '#000000'
            zorder = 5
 
        va = 'bottom' if y > 0 else 'top'
        ha = 'left' if x > 0 else 'right'
        xds = ds if x > 0 else -ds
        yds = ds if y > 0 else -ds
 
        txt = ax.annotate(name, xy=(x, y), alpha=text_alpha, color=text_color, xytext=(xds, yds), fontsize=label_size,
                          textcoords='offset points', ha=ha, va=va, clip_on=False, zorder=zorder)
        txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF80')])

    ax.axhline(nlp, **hline_kw)
    ax.axvline(nlp, **hline_kw)
    ax.axhline(-nlp, **hline_kw)
    ax.axvline(-nlp, **hline_kw)
    
    ax.set_xlabel(f'p-value {a1} vs {a2}')
    ax.set_ylabel(f'p-value {b1} vs {b2}')
    ax.set_ylim(-p_max, p_max)
    ax.set_xlim(-p_max, p_max)
        
    ax.scatter([], [], s=0.5*min_size, color='w', edgecolor='k', linewidth=0.5, label='low min. fold-change')
    ax.scatter([], [], s=5*min_size, color='w', edgecolor='k', linewidth=0.5, label='high min. fold-change')
    ax.scatter([], [], s=2*min_size, color=cmap(0.9), edgecolor='k', linewidth=0.5, label='hit +fold change')
    ax.scatter([], [], s=2*min_size, color=cmap(0.1), edgecolor='k', linewidth=0.5, label='hit -fold change')
    ax.scatter([], [], s=2*min_size, color=cmap(0.5), label='insignificant')
    if min_peps > 1:
       ax.scatter([], [], s=2*min_size, color=lcolor, edgecolor='none', alpha=0.5, label='low pep count')
    ax.scatter([], [], s=32, marker="$X^!$", color='#800080', edgecolor='none', label='comparison missing')
    ax.scatter([], [], s=32, marker="$X^0$", color='#400040', edgecolor='none', label='vs all zeros')
    ax.scatter([], [], s=32, marker="$X^1$", color='#400040', edgecolor='none', label='vs single')
              
    ax.legend(fontsize=6.2, loc='lower right')
    
    p = int(p_max)
    ticks = np.logspace(-p, 0, p+1, base=10.0)
    ticks = ticks[ticks > 2.0 ** (-p_max)]
    tpos = -np.log2(ticks)
    
    tpos = np.concatenate([tpos, -tpos[:-1][::-1]])
    ticks = np.concatenate([ticks, ticks[:-1][::-1]])
    ticks = _nice_tick_labels(ticks)
    
    ax.set_yticks(tpos)
    ax.set_yticklabels(ticks, fontsize=8)
    ax.set_xticks(tpos)
    ax.set_xticklabels(ticks, fontsize=8)
    
    texts = [f'{a1} > {a2}', f'{a2} > {a1}', f'{b1} > {b2}', f'{b2} > {b1}']
    max_len = max([len(x) for x in texts])
    fontsize = min(14, 280.0/max_len)
    
    # X    
    ax.text(1.0, -0.04, f'{a1} > {a2}', ha='right', va='top', color=cmap(1.0), transform=ax.transAxes, fontsize=fontsize, clip_on=False)
    ax.text(0.0, -0.04, f'{a2} > {a1}', ha='left',  va='top', color=cmap(0.0), transform=ax.transAxes, fontsize=fontsize, clip_on=False)
    
    # Y
    ax.text(-0.065, 1.0, f'{b1} > {b2}', ha='right', va='top', rotation=90, color=cmap(1.0), transform=ax.transAxes, fontsize=fontsize, clip_on=False)
    ax.text(-0.065, 0.0, f'{b2} > {b1}', ha='right',  va='bottom', rotation=90, color=cmap(0.0), transform=ax.transAxes, fontsize=fontsize, clip_on=False)
    
    xtlos = ax.get_xticklabels()
    ytlos = ax.get_yticklabels()
    for i, val in enumerate(tpos):
        if val > 0:
            xtlos[i].set_color(cmap(1.0))
            ytlos[i].set_color(cmap(1.0))
        if val < 0:
            xtlos[i].set_color(cmap(0.0))
            ytlos[i].set_color(cmap(0.0))
    
    plt.subplots_adjust(wspace=0.0, hspace=0.0, top=0.9, bottom=0.1, left=0.1, right=0.9)
    
    _watermark(fig)

    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()


def dual_comparison_plot(pdf, pair1, df1, pair2, df2, p_thresh, fc_thresh, min_peps, markers=None, is_pp=True,
                         colors=['#0080FF','#A0A000','#D00000'], max_size=160.0, min_size=8.0, label_size=5.0):

    # ? Add option to have these plots
    # Table output: fcs and pvals, max(pvalue) > thresh, max(p_value) > thresh (different colours), n_obs 1,2,3,4, n_pep, samce accession, labels, extra cols as volc.
    
    a1, a2 = pair1
    b1, b2 = pair2
    nlp = - np.log2(p_thresh)
    
    cmap = LinearSegmentedColormap.from_list('dual_plot', colors)
    grey = '#80808080'
    lcolor = '#808080'
    
    hline_kw = dict(ls='--', alpha=0.5, lw=0.5, color=lcolor)
    index1 = df1.index
    index2 = df2.index
    
    common = list(set(index1) & set(index2))
    not1 = list(set(index2) - set(index1))
    not2 = list(set(index1) - set(index2))
    
    pvalues1 = np.array(df1.loc[common, 'Tpvals'])
    fcvalues1 = np.array(df1.loc[common, 'grp1grp2_FC'])
    
    pvalues2 = np.array(df2.loc[common, 'Tpvals'])
    fcvalues2 = np.array(df2.loc[common, 'grp1grp2_FC'])
    
    names = np.array(df1.loc[common, 'labels'])
    
    has_zeros = (np.array(df1.loc[common, 'nobs_grp1']) == 0) | \
                (np.array(df1.loc[common, 'nobs_grp2']) == 0) | \
                (np.array(df2.loc[common, 'nobs_grp1']) == 0) | \
                (np.array(df2.loc[common, 'nobs_grp2']) == 0)
    
    low_qual =  (np.array(df1.loc[common, 'nobs_grp1']) < 2) | \
                (np.array(df1.loc[common, 'nobs_grp2']) < 2) | \
                (np.array(df2.loc[common, 'nobs_grp1']) < 2) | \
                (np.array(df2.loc[common, 'nobs_grp2']) < 2)
   
    if 'npeps' in df1:
        npeps = np.array(df1.loc[common, 'npeps'])
    else:
        npeps = np.full(len(names), min_peps)
    
    dual = np.ones(len(npeps), bool)
    
    if not1: # Not in first comparison
        n = len(not1)
        pvalues1 = np.concatenate([pvalues1, np.zeros(n)])
        fcvalues1 = np.concatenate([fcvalues1, np.zeros(n)])
        pvalues2 = np.concatenate([pvalues2, np.array(df2.loc[not1, 'Tpvals'])])
        fcvalues2 = np.concatenate([fcvalues2,np.array(df2.loc[not1, 'grp1grp2_FC'])])
        names = np.concatenate([names,np.array(df2.loc[not1, 'labels'])])
        has_zeros = np.concatenate([has_zeros, np.zeros(n, bool)])
        low_qual = np.concatenate([low_qual, np.zeros(n, bool)])
        dual = np.concatenate([dual, np.zeros(n, bool)])
        if 'npeps' in df1:
            npeps = np.concatenate([npeps, np.array(df2.loc[not1, 'npeps'])])
        else:
            npeps = np.concatenate([npeps, np.full(len(not1), min_peps)])

    if not2: # Not in second comparison
        n = len(not2)
        pvalues2 = np.concatenate([pvalues2, np.zeros(n)])
        fcvalues2 = np.concatenate([fcvalues2, np.zeros(n)])
        pvalues1 = np.concatenate([pvalues1, np.array(df1.loc[not2, 'Tpvals'])])
        fcvalues1 = np.concatenate([fcvalues1, np.array(df1.loc[not2, 'grp1grp2_FC'])])
        names = np.concatenate([names, np.array(df1.loc[not2, 'labels'])])
        has_zeros = np.concatenate([has_zeros, np.zeros(n, bool)])
        low_qual = np.concatenate([low_qual, np.zeros(n, bool)])
        dual = np.concatenate([dual, np.zeros(n, bool)])

        if 'npeps' in df1:
            npeps = np.concatenate([npeps, np.array(df1.loc[not2, 'npeps'])])
        else:
            npeps = np.concatenate([npeps, np.full(len(not2), min_peps)])

    p_best = np.maximum(pvalues1, pvalues2) # neg logs; bigger better
    f_best = np.where(np.abs(fcvalues1) > np.abs(fcvalues2), fcvalues1, fcvalues2)
    
    is_marker = np.array([x in set(markers) for x in names])
    is_hit = ((pvalues1 >= nlp) & (np.abs(fcvalues1) >= fc_thresh)) | ((pvalues2 >= nlp) & (np.abs(fcvalues2) >= fc_thresh)) 
    is_hit |= is_marker

    p_best /= p_best.max() # 0..1
    f_best /= np.abs(f_best).max() # 0..1
    
    weights = f_best # -1..1
    weights /= np.abs(weights).max()

    sizes = weights * weights
    sizes = min_size + (max_size-min_size) * sizes
    
    weights = (1.0 + weights) / 2.0 # 0..1
    
    xy_min = 0.0
    
    if is_pp:
        xy_max = math.ceil(max(pvalues1.max(), pvalues2.max()))
        x_values = pvalues1
        y_values = pvalues2
        x_values[fcvalues1 < 0.0] *= -1.0
        y_values[fcvalues2 < 0.0] *= -1.0

    else:
        xy_max = math.ceil(max(np.abs(fcvalues1).max(), np.abs(fcvalues2).max()))
        x_values = fcvalues1
        y_values = fcvalues2
    
    s0 = (npeps < min_peps) & dual 
    s1 = ~is_hit & ~s0 & dual & ~is_marker
    s2 = is_hit & ~s0 & dual & ~is_marker
    s3 = is_hit & ~s0 & ~dual & ~is_marker
    
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_title(f'[{a1} vs {a2}] vs [{b1} vs {b2}]', fontsize=14)
    
    ax.scatter(x_values[s0], y_values[s0], color=lcolor, alpha=0.3, clip_on=False, s=sizes[s0], zorder=1)
    ax.scatter(x_values[s1], y_values[s1], c=[cmap(w) for w in weights[s1]], alpha=0.3, clip_on=False, s=sizes[s1], zorder=2)
    ax.scatter(x_values[s2], y_values[s2], c=[cmap(w) for w in weights[s2]], alpha=1.0, edgecolor='k', linewidth=0.5, clip_on=False, s=sizes[s2], zorder=4)
    ax.scatter(x_values[s3], y_values[s3], c=[cmap(w) for w in weights[s3]], alpha=1.0, edgecolor='#800080', linewidth=0.5, clip_on=False, s=sizes[s3], zorder=3)
    ax.scatter(x_values[is_marker], y_values[is_marker], c='k', alpha=1.0, marker='*', clip_on=False, s=1.5*sizes[is_marker], zorder=4)
    
    hit_idx = np.nonzero(is_hit)[0]
    for i in hit_idx:
        ds = 0.5 * np.sqrt(sizes[i])
        name = names[i]
        text_alpha = 1.0
        low_pep = s0[i]
        x = x_values[i]
        y = y_values[i]
 
        if not isinstance(name, str):
            name = '?'
 
        if low_pep:
            text_color = lcolor
            text_alpha = 0.5
            zorder = 1
        
        if not dual[i]:
            name = name + '$^{!}$'
            text_color = '#800080'
            text_alpha = 0.5
            zorder = 2
          
        elif has_zeros[i]:
            name = name + '$^{0}$'
            text_color = '#400040'
            zorder = 3
 
        elif low_qual[i]:
            name = name + '$^{1}$'
            text_color = '#400040'
            zorder = 3
 
        else:
            text_color = '#000000'
            zorder = 5
 
        va = 'bottom' if y > 0 else 'top'
        ha = 'left' if x > 0 else 'right'
        xds = ds if x > 0 else -ds
        yds = ds if y > 0 else -ds
 
        txt = ax.annotate(name, xy=(x, y), alpha=text_alpha, color=text_color, xytext=(xds, yds), fontsize=label_size,
                          textcoords='offset points', ha=ha, va=va, clip_on=False, zorder=zorder)
        txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF80')])

    line_val = nlp if is_pp else fc_thresh
    ax.axhline(line_val, **hline_kw)
    ax.axvline(line_val, **hline_kw)
    ax.axhline(-line_val, **hline_kw)
    ax.axvline(-line_val, **hline_kw)
    
    sz_tag = 'fold change' if is_pp else 'p-value'
    ax_tag = 'fold change' if (not is_pp) else 'p-value'

    ax.set_xlabel(f'{ax_tag} {a1} vs {a2}')
    ax.set_ylabel(f'{ax_tag} {b1} vs {b2}')
    ax.set_ylim(-xy_max, xy_max)
    ax.set_xlim(-xy_max, xy_max)
    
    ax.scatter([], [], s=0.5*min_size, color='w', edgecolor='k', linewidth=0.5, label=f'low min. {sz_tag}')
    ax.scatter([], [], s=5*min_size, color='w', edgecolor='k', linewidth=0.5, label=f'high min. {sz_tag}')
    ax.scatter([], [], s=2*min_size, color=cmap(0.9), edgecolor='k', linewidth=0.5, label='hit +fold change')
    ax.scatter([], [], s=2*min_size, color=cmap(0.1), edgecolor='k', linewidth=0.5, label='hit -fold change')
    ax.scatter([], [], s=2*min_size, color=cmap(0.5), label='insignificant')
    if min_peps > 1:
       ax.scatter([], [], s=2*min_size, color=lcolor, edgecolor='none', alpha=0.5, label='low pep count')
    if markers:
       ax.scatter([], [], s=2*min_size, color='k', marker='*', label='marker')
    ax.scatter([], [], s=32, marker="$X^!$", color='#800080', edgecolor='none', label='comparison missing')
    ax.scatter([], [], s=32, marker="$X^0$", color='#400040', edgecolor='none', label='vs all zeros')
    ax.scatter([], [], s=32, marker="$X^1$", color='#400040', edgecolor='none', label='vs single')
              
    ax.legend(fontsize=6.2, loc='lower right')
    
    p = int(xy_max)
    
    if is_pp:
        ticks = np.logspace(-p, 0, p+1, base=10.0)
        ticks = ticks[ticks > 2.0 ** (-xy_max)]
        tpos = -np.log2(ticks)
    else:
        ticks = np.logspace(p, 0, p+1, base=2.0)
        ticks = ticks[ticks < 2.0 ** (xy_max)]
        tpos = np.log2(ticks)

    tpos = np.concatenate([tpos, -tpos[:-1][::-1]])
    ticks = np.concatenate([ticks, ticks[:-1][::-1]])
    ticks = _nice_tick_labels(ticks)
    
    ax.set_yticks(tpos)
    ax.set_yticklabels(ticks, fontsize=8)
    ax.set_xticks(tpos)
    ax.set_xticklabels(ticks, fontsize=8)
    
    texts = [f'{a1} > {a2}', f'{a2} > {a1}', f'{b1} > {b2}', f'{b2} > {b1}']
    max_len = max([len(x) for x in texts])
    fontsize = min(14, 280.0/max_len)
    
    # X    
    ax.text(1.0, -0.04, f'{a1} > {a2}', ha='right', va='top', color=cmap(1.0), transform=ax.transAxes, fontsize=fontsize, clip_on=False)
    ax.text(0.0, -0.04, f'{a2} > {a1}', ha='left',  va='top', color=cmap(0.0), transform=ax.transAxes, fontsize=fontsize, clip_on=False)
    
    # Y
    ax.text(-0.065, 1.0, f'{b1} > {b2}', ha='right', va='top', rotation=90, color=cmap(1.0), transform=ax.transAxes, fontsize=fontsize, clip_on=False)
    ax.text(-0.065, 0.0, f'{b2} > {b1}', ha='right',  va='bottom', rotation=90, color=cmap(0.0), transform=ax.transAxes, fontsize=fontsize, clip_on=False)
    
    xtlos = ax.get_xticklabels()
    ytlos = ax.get_yticklabels()
    for i, val in enumerate(tpos):
        if val > 0:
            xtlos[i].set_color(cmap(1.0))
            ytlos[i].set_color(cmap(1.0))
        if val < 0:
            xtlos[i].set_color(cmap(0.0))
            ytlos[i].set_color(cmap(0.0))
    
    plt.subplots_adjust(wspace=0.0, hspace=0.0, top=0.9, bottom=0.1, left=0.1, right=0.9)
    
    _watermark(fig)

    if pdf:
        pdf.savefig(dpi=DPI)
    else:
        plt.show()
    
    plt.close()


def volcano(pdf, pair_name, df, FClim, pthresh, min_peps, colors=['#0080FF','#A0A000','#D00000'],
            quant_scale=False, split_x=True, hq_only=False, hit_labels=True, markers=None, lw=0.25,
            ls='--', lcolor='#808080', max_size=160.0, min_size=8.0, label_size=5.0, show_low_pep=True, marker_text_col=None):
    
    cmap = LinearSegmentedColormap.from_list('volc', colors[::-1])
    group1, group2 = pair_name.split(':::')
    
    if not markers:
        markers = []
        
    pthresh_orig = pthresh
    pthresh = -np.log2(pthresh)
    
    marker_texts = df[marker_text_col] if marker_text_col else None
    
    high_qual = (np.array(df['nobs_grp1']) > 1) & (np.array(df['nobs_grp2']) > 1)
    has_zeros = (np.array(df['nobs_grp1']) == 0) | (np.array(df['nobs_grp2']) == 0)
    intensity = np.maximum(np.array(df['zmean_grp1']), np.array(df['zmean_grp2']))
    names = np.array(df['labels'])
    

    if 'npeps' in df:
        npeps = np.array(df['npeps'])
    else:
        npeps = np.full(len(names), min_peps)

    pvalues = np.array(df['Tpvals'])
    fcvalues = np.array(df['grp1grp2_FC'])
    
    if hq_only:
        pvalues = pvalues[high_qual]
        fcvalues = fcvalues[high_qual]
        names = names[high_qual]
        high_qual = high_qual[high_qual]
        has_zeros = high_qual[high_qual]
 
    color_pos, color_low, color_neg = colors
    
    if split_x:
         xlims, neg_split, pos_split = util.get_splitxlist(df, abs(split_x))
    
    else:
         xlims = util.get_xlist(df)
         neg_split, pos_split = None, None
    
    y_label = 'p-value'
    p_label = f'{pthresh_orig:.3g}'
    p_label_kw = dict(xytext=(2,2), xycoords='data', textcoords='offset points',
                      color='#808080', fontsize=5, ha='left', va='bottom')
    
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
    
    if marker_texts is None:
        title = pair_name.replace(':::',' vs ')
    else:
        title = pair_name.replace(':::',' vs ')
        title = 'Max. observed peptides for ' + title
    
    ax1.set_title(title, fontsize=18)
    ax1.set_xlabel('log2 fold change')
    
    pnorm = np.array(pvalues)
    pnorm /= pnorm.max() # 0..1
    
    fcnorm = np.array(fcvalues)
    fcnorm /= np.abs(fcnorm).max() # 0..1
    
    weights = pnorm * fcnorm # -1..1
    weights /= np.abs(weights).max()
    weights = (1.0 + weights) / 2.0 # 0..1
    
    if quant_scale:
        intensity = np.argsort(intensity).astype(float)
    else:
       intensity -= intensity.min()
       
    intensity /= intensity.max()

    #plt.hist(intensity, bins=100)
    #plt.show()

    intensity **= 2.0
    
    sizes = min_size + (max_size-min_size) * intensity # min_size .. max_size
    
    scatter_kw = dict(alpha=0.3, edgecolors='none', clip_on=False)
    hit_kw = dict(alpha=1.0, clip_on=False, edgecolor='k', linewidth=0.5)
    hit_txt_kw = dict(alpha=1.0, clip_on=False, va='center', ha='center', fontsize=label_size)
   
    # middle bounds
    m0 = xmin
    m1 = xmax
    
    unimportant = (pvalues < pthresh) | (np.abs(fcvalues) < FClim) | (npeps < min_peps)
    
    if ax0:
        selection = (fcvalues < neg_split[0]) & high_qual & unimportant
        
        if np.any(selection):
            colors = [cmap(w) for w in weights[selection]]
            ax0.scatter(fcvalues[selection], pvalues[selection], s=sizes[selection], c=colors,**scatter_kw)
        
        m0 = neg_split[1]
        
    if ax2:
        selection = (fcvalues > pos_split[1])  & high_qual & unimportant
        if np.any(selection):
            colors = [cmap(w) for w in weights[selection]]
            ax2.scatter(fcvalues[selection], pvalues[selection], s=sizes[selection], c=colors, **scatter_kw)
        
        m1 = pos_split[0]
    
    selection = (fcvalues >= m0) & (fcvalues <= m1) & high_qual & unimportant
    
    if np.any(selection):
        colors = [cmap(w) for w in weights[selection]]
        ax1.scatter(fcvalues[selection], pvalues[selection], s=sizes[selection], c=colors, **scatter_kw)

    marker_text = []
    
    for i, name in enumerate(names):
        size = sizes[i]
        low_qual = not high_qual[i]
        color = cmap(weights[i])
        ds = 0.5 * np.sqrt(size)
        x = fcvalues[i]
        y = pvalues[i]
        
        if ax0 and (x <= neg_split[0]):
            ax = ax0
        elif ax2 and (x >= pos_split[1]):
            ax = ax2
        else:
            ax = ax1  
 
        if name in markers:
            marker_text.append((x, y, ds, name))
            ax.scatter(x, y, s=1.5*size, color='k', marker='*', edgecolor='k', clip_on=False, zorder=11)
            continue
        
        low_pep = npeps[i] < min_peps
        
        if low_pep and not show_low_pep:
            continue
        
        if hq_only and low_qual:
            continue
        
        if marker_texts is None:
            text_alpha = 1.0
            
            if not isinstance(name, str):
                name = '?'
            
            if low_pep:
                text_color = lcolor
                text_alpha = 0.5
               
            elif has_zeros[i]:
                name = name + '$^{0}$'
                text_color = '#400040'
 
            elif low_qual:
                name = name + '$^{1}$'
                text_color = '#400040'
 
            else:
                text_color = '#000000'
        else:
            text_color = '#00000080'
        
        pos_hit = x >= FClim and y >= pthresh
        neg_hit = x <= -FClim and y >= pthresh
        
        
        if pos_hit or neg_hit:
            if marker_texts is None:
                if low_pep:
                    z_order = 9
                    ax.scatter(x, y, s=size, color=lcolor, edgecolor='none', alpha=0.5, zorder=z_order)
                else:
                    z_order = 11
                    ax.scatter(x, y, s=size, color=color, zorder=size, **hit_kw)
            else:
                z_order = 11
                ax.text(x, y, marker_texts[i], color='k', zorder=z_order, **hit_txt_kw)
                text_color = color
                
            if hit_labels:
                if pos_hit:
                    xytext=(ds, ds)
                    ha='left'
                
                else:
                    xytext=(-ds, ds)
                    ha='right'
                
                for a in (ax0, ax1, ax2):
                    if not a:
                        continue
 
                    txt = a.annotate(name, xy=(x, y), alpha=text_alpha, color=text_color, xytext=xytext, fontsize=label_size,
                                     textcoords='offset points', ha=ha, clip_on=False, zorder=z_order+1)
                    txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF80')])
             
            continue
             
        if low_qual:
            ax.scatter(x, y, s=min_size, alpha=0.5, edgecolor='none', color=color, zorder=9) # color=color_low           
               
    for x, y, ds, name in marker_text:
         txt = ax.annotate(name, xy=(x, y), xytext=(ds, ds), textcoords='offset points',
                          fontsize=label_size, clip_on=False)
         txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='#FFFFFF80')])
    
    ax = ax0 if ax0 else ax1
    ax.text(0.05, -0.05, f'Greater in {group2}', color=color_neg, ha='left', va='top', transform=ax.transAxes, fontsize=10, clip_on=False)
    
    ax.set_ylim(0.0)
    pmax = math.ceil(np.nanmax(pvalues))
    p = int(pmax) 
    y_ticks = np.logspace(-p, 0, p+1, base=10.0)
    y_ticks = y_ticks[y_ticks > 2.0 ** (-pmax)]
    
    ax.yaxis.tick_left()
    ax.set_yticks(-np.log2(y_ticks))
    ax.set_yticklabels(_nice_tick_labels(y_ticks))
    
    ax = ax2 if ax2 else ax1
    ax.text(0.95, -0.05, f'Greater in {group1}', color=color_pos, ha='right', va='top', transform=ax.transAxes, fontsize=10, clip_on=False)
     
    ax.scatter([], [], s=0.5*min_size, color='w', edgecolor='k', linewidth=0.5, label='low abundance')
    ax.scatter([], [], s=5*min_size, color='w', edgecolor='k', linewidth=0.5, label='high abundance')
    ax.scatter([], [], s=2*min_size, color=cmap(0.9), edgecolor='k', linewidth=0.5, label='+ve hit')
    ax.scatter([], [], s=2*min_size, color=cmap(0.1), edgecolor='k', linewidth=0.5, label='-ve hit')
    ax.scatter([], [], s=2*min_size, color=cmap(0.5), label='insignificant')
    if min_peps > 1:
       ax.scatter([], [], s=2*min_size, color=lcolor, edgecolor='none', alpha=0.5, label='low pep count')
    if markers:
       ax.scatter([], [], s=2*min_size, color='k', marker='*', label='marker')
    ax.scatter([], [], s=32, marker="$X^0$", color='#400040', edgecolor='none', label='vs all zeros')
    ax.scatter([], [], s=32, marker="$X^1$", color='#400040', edgecolor='none', label='vs single')
              
    ax.legend(fontsize=6.2, loc='lower right')
    
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
    

"""
gridsize = 100
colors=['#00FFFF','#0000FF','#000000','#A00000','#FFFF00']
show_pairs = defaultdict(list)
for g1, g2 in pairs:
    bg1 = bg_dict.get(g1, 'None')
    bg2 = bg_dict.get(g2, 'None')
    
    if bg1 != bg2:
        key = tuple(sorted([bg1, bg2]))
        show_pairs[key].append((g1, g2))

n_rows = len(show_pairs)
n_cols = 2 # before, after
fig, axarr = plt.subplots(n_rows, n_cols)
cmap = LinearSegmentedColormap.from_list(colors)

for row, key in enumerate(show_pairs):
    bg1, bg2 = key
    x_vals1 = []
    y_vals1 = []
    x_vals2 = []
    y_vals2 = []
    for g1, g2 in show_pairs[key]:
        x_vals1.append(df[g1])
        y_vals1.append(df[g2])
        x_vals2.append(df['znorm_'+g1])
        y_vals2.append(df['znorm_'+g2])
      
    ax1 = axarr[row, 0]
    ax2 = axarr[row, 1]
    
    ax1.hexbin(x_vals1, y_vals1, cmap=cmap, gridsize=gridsize, mincnt=1, edgecolors='none')
    ax2.hexbin(x_vals2, y_vals2, cmap=cmap, gridsize=gridsize, mincnt=1, edgecolors='none')
    
    ax1.set_xlabel('{bg1} $log_2$ count')
    ax1.set_ylabel('{bg2} $log_2$ count')
    ax2.set_xlabel('{bg1} Adj. $log_2$ count')
    ax2.set_ylabel('{bg2} Adj. $log_2$ count')
    
    for row == 0:
      ax1.set_title('Original comparisons (aggregated)')
      ax2.set_title('Adjusted comparisons')
    
"""
