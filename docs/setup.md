# Setup  

For this tutorial, there are three components to set up. Depending on what you'd like to do, you may need one, two, or all of these components.  


## 1. Make a copy of this tutorial [required for all]  

There are two ways to do this:
* [Recommended] If you're familiar with git, clone this repository either via the web interface, a GUI such as gitkraken, or the command line:  
```bash  
git clone https://github.com/nextstrain/ncov.git
```

* [Alternative] If you're not familiar with git, you can also download a copy of these files by navigating to [https://github.com/nextstrain/ncov](https://github.com/nextstrain/ncov) and clicking on the green 'Clone or download' button.


## 2. Install augur [required for running analysis]  

You can find instructions for installing augur [here](https://nextstrain.org/docs/getting-started/introduction).

Note there are many options for installing augur, including local installations from `conda`, `pip`, and from source, and container installation - take a minute to look at which may be best for you.

If you are running Windows, we recommend working through Windows Subsystem for Linux (WSL) - you can find out how to set this up [here](https://nextstrain.org/docs/getting-started/windows-help).


## 3. Install auspice [required only for local visualization; web-based, no-install options available]

You can find instructions for installing auspice [here](https://nextstrain.github.io/auspice/introduction/install).

---

## Advanced reading: considerations for keeping a 'Location Build' up-to-date

_Note: we'll walk through what each of the referenced files does shortly_

If you are aiming to create a public health build for a state, division, or area of interest, you likely want to keep your analysis up-to-date easily.
If your run contains contextual subsampling (sequences from outside of your focal area), you should first ensure that you [regularly download the latest sequences as input](data-prep.md), then re-run the build.
This way, you always have a build that reflects the most recent SARS-CoV-2 information.

You should also aim to keep the `ncov` repository updated.
If you've clone the repository from Github, this is done by running `git pull`.
This downloads any changes that we have made to the repository to your own computer.
In particular, we add [new colors and latitute & longitude information](customizing-analysis.md) regularly - these should match the new sequences you download, so that you don't need to add this information yourself.

If you don't need to share the contents of [`my_config`](orientation-files.md) with anyone, then you can leave this in the `./my_config/` folder.
It won't be changed when you `git pull` for the latest information.

However, if you want to share your build config, you'll need to adopt one of the following solutions.
First, you can 'fork' the entire `ncov` repository, which means you have your own copy of the repository.
You can then add your config files to the repository and anyone else can download them as part of your 'fork' of the repository.
Note that if you do this, you should ensure you `pull` regularly from the original `ncov` repository to keep it up-to-date.

Alternatively, you can create a new repository to hold your `my_config` files, outside of the `ncov` repository.
You can then share this repository with others, and it's straightforward to keep `ncov` up to date, as you don't change it at all.
If doing this, it can be easiest to create a `my_config` folder and imitate the structure found in the `my_config` folder within `ncov`, but this isn't required.
Note that to run the build you'll need still run the `snakemake` command from within the `ncov` repository, but specify that the build you want is outside that folder.

For the [`south-usa-sarscov2`](https://github.com/emmahodcroft/south-usa-sarscov2/) example, you can see the `south-central` build set up in a `profiles` folder.
To run this, one would call the following from within `ncov`:

```bash
ncov$ snakemake --profile ../south-usa-sarscov2/profiles/south-central/
```


## [Previous Section: Preparing your data](./docs/data-prep.md)
## [Next Section: Orientation: analysis workflow](./docs/orientation-workflow.md)
