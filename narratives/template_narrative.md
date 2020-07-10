---
title: Genomic analysis of COVID-19 spread.
authors:
  - Name 1
  - Name 2

authorLinks:
  - https://author1.org
  - https://author2.io
affiliations: "Fred Hutch, Seattle, USA; Biozentrum, Basel, Switzerland; CZI, CA, USA"

license: "CC-BY"  
licenseLink: "https://creativecommons.org/licenses/by/4.0/"
dataset: "https://nextstrain.org/ncov/global?legend=closed" # must be accessible to the auspice server running the narrative

abstract: "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
---

<!-- Comment tags like these are not rendered, they're just helpful for you -->
<!-- Known 'gotcha' bug: ensure that links always end in a 'letter' (a period counts). If some kind of text doesn't follow them, it breaks the slide. -->


<!-- ############ SLIDE BREAK ############# -->
<!-- SLIDE 1 -->
<!--  Each slide MUST start with a link to a specific view of the dataset (must match the `dataset` specified above) -->
# [SLIDE 1 TITLE](https://nextstrain.org/ncov/global?c=country)

<!-- This is left-side text -->
Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.

<!-- There is NO right-side text on this slide -->


<!-- ############ SLIDE BREAK ############# -->
<!-- SLIDE 2 -->
# [SLIDE 2 TITLE](https://nextstrain.org/ncov/global?c=region)

<!-- This is the left-side text -->

[Including a link as an example; always end the line with a period or other 'letter' character](google.com)!
Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# Example using markdown for the right side  
We can also replace the right side view with whatever markdown contents we choose, including links, images, etc.
```
