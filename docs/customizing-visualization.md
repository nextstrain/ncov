# Customizing your Auspice visualization  

Just as we can specify a build-specific analysis options in the `config.yaml` file, we can also specify build-specific visualization options in this directory. An example of this can be seen in the `./my_config/example_advanced_customization/` directory.

Looking at the `config.yaml` file, the last few lines are:  
```yaml
files:
  colors: "my_config/example_advanced_customization/colors.tsv"
  auspice_config: "my_config/example_advanced_customization/auspice_config_swiss.json"
```

Let's look at what kinds of customization options we can use these for.


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
  colors: "my_config/<name>/colors.tsv"
```
to your `config.yaml` file.


## Changing the dataset description

The dataset description, which appears below the visualizations, is specified in `default_config/description.md`

## Adding custom metadata fields to color by   
1. Add a [valid metadata column](data-prep.md) to your `metadata.tsv`  
2. Using our example as a template, create an `auspice_config.json` file to your analysis directory
`sarscov2-tutorial$ cp ./my_config/example_advanced_customization/auspice_config_swiss.json ./my_config/<name>/`    
3. Add an entry to the `colorings` block of this JSON. You can see an example of this in `my_config/example_advanced_customization/auspice_config_swiss.json` file:  
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

## [Previous Section: Customizing your analysis](customizing-analysis.md)
## [Next Section: Options for visualizing and sharing results](sharing.md)
