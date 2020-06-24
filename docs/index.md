# Getting started

This template and tutorial will walk you through the process of running a basic phylogenetic analysis on SARS-CoV-2 data.
We've created these resources with the goal of enabling Departments of Public Health to start using Nextstrain to understand their SARS-CoV-2 genomic data within 1-2 hours.

To get started, [head over to the first section on preparing your data](data-prep.md)!

## Help  

If something in this tutorial is broken or unclear, please [open an issue](https://github.com/nextstrain/ncov/issues/new/choose) so we can improve it for everyone.  

If you have a specific question, post a note over at the [discussion board](https://discussion.nextstrain.org/) -- we're happy to help!


## Quickstart    

If you'd prefer, you can also start by running a basic analysis on the provided example data and/or visualizing the output with the [auspice.us](auspice.us) drag-and-drop viewer. If you get stuck at any point, you can find more detailed instructions in the full tutorial outlined above.

We also recommend [this 1-hour video overview](https://youtu.be/m4_F2tG58Pc) by Heather Blankenship on how to deploy Nextstrain for a Public Health lab.

#### 1. Clone this repository  
```
git clone https://github.com/nextstrain/ncov.git
```

#### 2. Install augur for bioinformatics
```
python3 -m pip install nextstrain-augur
```

MacOS:
```
brew install mafft iqtree raxml fasttree vcftools
```

Ubuntu/Debian:  
```
sudo apt install mafft iqtree raxml fasttree vcftools
```

#### 3. Run basic analysis on example data  
```
ncov$ snakemake --profile ./my_config/example
```


#### 4. Visualize your results (or our example output)  
Go to `https://auspice.us` in your browser.
Drag and drop `./auspice/sarscov2.json` (or any other JSON in this directory) anywhere on the screen.

Voila!
