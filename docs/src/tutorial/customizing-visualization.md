# Customizing your Auspice visualization  

Just as we can specify a build-specific analysis options in the `builds.yaml` file, we can also specify build-specific visualization options in this directory.

Looking at the `builds.yaml` file, the last few lines are:  
```yaml
files:
  auspice_config: "my_profiles/example/my_auspice_config.json"
```

This points to a JSON file that parameterizes the output files used for visualizion with Auspice.
Let's look at what kinds of customization options we can use this for.


## Custom color schemes  

If you'd like to specify a custom color scale, you can add a `colors.tsv` file, where each line is a tab-delimited list of a metadata column name; a metadata value; and a corresponding hex code.   

The first few lines of the example file look like this:  
```
country	Russia	#5E1D9D
country	Serbia	#4D22AD
country	Europe	#4530BB
...
```

Make sure to also add
```yaml
files:
  colors: "my_profiles/<name>/colors.tsv"
```
to your `builds.yaml` file.


## Changing the dataset description

The dataset description, which appears below the visualizations, is read from a file which is specified in `builds.yaml`.  Per-build description can be set by specifying them in the build.

```yaml
builds:
  north-america: # name of the build; this can be anything
    description: my_profiles/example/north-america-description.md
```

If that is not provided, then a per-run description is used, also specified in `builds.yaml`:

```yaml
files:
  description: my_profiles/example/my_description.md
```

## Adding custom metadata fields to color by   
1. Add a [valid metadata column](data-prep.md) to your `metadata.tsv`  
2. Open `my_profiles/<name>/auspice_config.json`  
3. Add an entry to the `colorings` block of this JSON:

```json
...
"colorings": [
  {
    "key": "location",
    "title": "Location",
    "type": "categorical"
  },
  {
    "key": "metadata_column_name",
    "title": "Display name for interface",
    "type": "categorical" or "continuous"
  }
...
]
...
```

## Choosing defaults  
You can specify the default view in the `display_defaults` block of an `auspice_config.json` file (see above)
```json
...
"display_defaults": {
  "color_by": "division",
  "distance_measure": "num_date",
  "geo_resolution": "division",
  "map_triplicate": true,
  "branch_label": "none"
},
...
```

## Choosing panels to display  

Similarly, you can choose which panels to enable in the `panels` block:  
```json
...
"panels": [
  "tree",
  "map",
  "entropy"
]
...
```

