# LAVA: a minimal-imputational method for dealing with missing values in experimental with low replicate numbers

## Introduction  

Welcome to the Lava repository. Lava works on Linux, Windows and Mac (though it's not been tested on a Mac for a long time).  
This software has been built for the analysis of exploratory proteomics experiments that have few replicates <strong>(minimum 3)</strong>. Missing data values can be much more problematic in proteomics compared to sequencing data and if a pilot study only contains a few replicates, thye way the missing values are handled can lead to dramatically different conclusions. Some off-the-shelf differential abundance software can miss out proteins where that value is missing in one group of replicates. Or, more worryingly, impute that value from the other group. Imagine that this protein had been knocked out in one group, so you knew it should be absent from all those replicates. That could be the most important conslusion, so you don't want the differences flattened out by imputation. You might want to set missing values to the lowest detectable value of the whole dataset. That's a fair idea, if the software allows it, but is this going to affect the variance? The number of missing values and the variance of the unimputed values are important features of any dataset, so need to be preserved.  

Lava tried to preserve the maximum number of data points and their variance within any given dataset. It does this by accounting for replicate number in each treatment group and applying these rules:  

<li> Half or more values present in both groups: measurement accepted  
<li> All values absent in one group and over half present in corresponding group: measurement accepted  
<li> Single value present in one group and all values present in corresponding group: measurement <strong>optionally</strong> accepted
<li> Fewer than half, but more than zero, values present in both groups: measurement rejected  

<br>  

These rules triage the data, keeping only those where there's enough confidence to estimate the true value of the measurement. For any measurement that's been accepted, as long as it has at least two values in each group, we can perform a T-test for unequal group sizes using ```scipy.stats.ttest_ind_from_stats``` , thus preserving the original variance. When all values are absent in one group and mostly present in the other, we consider this good evidence that the protein really is below the limit of detection in one group. We set this value to just below the lowest reported, and set its variance to the mean of the corresponding group. These datapoints are clearly highlighted in output files, allowing users to see (or ignore) their impact. How abuot when a single value is present in one group and most are present in the other group? This touches on the balance between having quality results vs. seeing what might be happening in a preliminary study, so we let users decide with the -hq (high quality) argument. Omitting this argument means that, in the group where a single value was present, this value is taken as the mean but the data is penanlised with a variance at the 80th centile for the dataset. Naturally, these datapoints are also highlighted in the outputs.  

<br>  

Lava also includes other features that, as a proteomics researcher, I've found useful. Results can be filtered to remove proteins with fewer than a set number of peptides. Labels can be to Uniprot ID or Gene name, marker proteins highlighted on plots, p-value and fold-change thresholds altered.  

### Fold-Change Fold-Change plots and P-value P-value plots  
These plots are unique to Lava. They compare either the Fold-Change or P-value output of one pair against another. For example, if you want to be looking at the fold change of a protein between mutant1 and mutant2, but each mutant must also be compared to a control, you can plot fold-change (control-mutant1) vs foldchange (control-mutant2) and get all fold-change comparisons in one plot. Values for the unplotted data (P-values, if plotting a Fold-Change plot) are included in point colour, so there's no loss of information.  

We have tried to make aesthetically pleasing plots that are ready for publication but changes can be made in the plot.py file.  


## How to use Lava (scroll down for troubleshooting)

Assuming you've got git or git-bash for Windows etc. installed, and you are already at https://github.com/tempeparsons/Lava, navigate via the terminal to your relevant directory, then paste this into the terminal:  

```git clone https://github.com/tempeparsons/Lava.git ```  
 
Once the Lava repository has been cloned, you should see a directory structure like so:  

```
Lava    
├── scripts    
│   │   
│   └── spec_count_psm.py  
├── test_data  
│   │ 
│   ├── background_table.csv  
│   ├── experiment_table_dianndata.csv  
│   ├── experimental_table_PD3data.csv  
│   ├── VolcCLI_DIANN.txt    
│   └── VolcCLI_PD3.txt 
│
├── lava_utils.py  
├── lava_plots.py  
├── lava.py 
└── lava_environment.yml  
```  

The first step is to create the environment:  
(Lava only requires any recent version of numpy, pandas, matplotlib, scipy and scikit-learn, so you can build your preferred environment easily enough)  
```
conda create --file lava_environment.yml  
conda activate lava_env   
```  

You are now ready to start using the Lava software. 

Typing ```python3 lava.py -h``` lists and explains all the options for running Lava.   

Now try typing in the following command:  

```
>python3 lava.py test_data/VolcCLI_DIANN.txt -e test_data/experimental_table_dianndata.csv -i Protein.Group -g test_results_diann.pdf -o test_results_diann.csv
```

In the directory containing lava.py, a pdf and csv file should have been written. 

The .pfg graphical output contains QC plots and volcano plots showing an everything-by-everything comparison of the groups from VolcCLI_DIANN.txt that were included in experimental_table_2.csv. If you open up experimental_table_2.csv in the test_data folder, notice that not all possible groups need be included, and others can be commented outThe xy coordinates of each datapoint can be in the numerical output .csv file. Note the -i command specifying the non-default protein accession column name.  

Other things to notice in the graphical output are that your input arguments are stated at the top of each graphical output, so there can be no confusion over parameter choice. Also note the legends in the volcano plots. Some points are annotated as 'vs all zeros' or 'vs single'. This relates to how Lava manages to preserve high and medium quality data points without resorting to imputation, as described earlier. Data points are preserved at the user's discretion; try including -hq or --hq-only as an argument and see how the'vs single' datapoints are no longer plotted. 

Here's a slightly more complex command that demonstrates further utilities within Lava:

```
python3 lava.py test_data/VolcCLI_DIANN.txt -e test_data/experimental_table_dianndata.csv -i Protein.Group -g test_results_diann_bckg_mrks.pdf -o test_results_diann_bckg_mrks.csv -b test_data/background_table.csv -m 'P21333' 'Q8NF91' 'P55196' 'Q14108'
```
This command uses 'test_data/background_table.csv' to identify and subtract background control data from the experimental data. It also shows how select markers, passed to the command, can be highlighted as black stars on the plots. 

Next, try this:

```
python3 lava.py test_data/VolcCLI_DIANN.txt -e test_data/experimental_table_dianndata.csv -i Protein.Group -ff -pp -g test_results_diann_ffpp.pdf -o test_results_diann_ffpp.csv
```  
This demonstrates the -pp and -ff commands. They make plots of p-value vs p-value or fold-change vs fold-change, allowing you make comparisons between the output of a pair of pairs.  

Lava also allows the user to filter their data by the number of peptides associated with a protein. For the highest-confidence results only, the minimum peptide setting can be set to 2 (or higher). Of course, this only works if you have a peptide numbers column in your data. Depending on which software you've searched your mass spectrometry data with, you might need to merge this column in from a peptide table. Note how in the graphical output you can still see proteins identified by less than the minimum peptide threshhold, but they're greyed-out. Note also that Proteome Discoverer may put special characteris in columns names this peptides column name has been wrapped in quotation marks.  

```
python3 lava.py test_data/VolcCLI_PD3.txt -e test_data/experimental_table_PD3data.csv -n 2 -nc '# Peptides' -g test_results_PD3.pdf -o test_results_PD3.csv
```


TMT-labelled data: This works fine with Lava, but you may find data points sizes come out looking very similar if you're plotting ratios of isobaric tags . The ```-qs, --quantile-scaling``` argument has been designed for this; it sets the datapoint sizing scale relative to your dataset, rather than just setting point size to mean abudnance value.  

This README intriduced you to the basic concepts of Lava. Experiment with the other arguments and please get in touch if anything isnt' working.  

## Troubleshooting the install  

You can check the install with by running ```python check_install.py```.  

The conda install has been tested on Linux and Windows (assuming Anaconda or Miniconda installed), but if used within an institute, firewalls may block it. If you're using windows with miniconda and getting access error messages, try calling the environment's einterpreter directly. This will probably looking something like:

```C:\Users\<your_user_name>\AppData\Local\miniconda3\envs\lava_env\python.exe lava.py <rest of the command> ```


Alternatively, use the requirements.txt file from Windows Powershell:  
```
python -m venv lava_env
.\lava_env\Scripts\python.exe -m pip install -r lava_requirements.txt
.\lava_env\Scripts\python.exe lava.py --help
```  

You can cehck the install by running:  
```
.\lava_env\Scripts\python.exe check_install.py
.\lava_env\Scripts\python.exe lava.py --help
```


