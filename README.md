### Intro
BacDive is the largest aggregate database for bacterial information. From here, we’ve built a visualization tool that allows the user to see resistance, sensitivity, production, and utilization trends between strains, species, genera, and more. The Bacdive Streamlit can be used to generate a great number of graphs and tables. The user by default gets information on metabolite resistances, production, utilization, and sensation. The user can further select whether they want to view strains of bacteria species, or exclude them. Furthermore, instead of only showing positive results, negative results can also be chosen. It can detail the shared antibiotic resistances, production, utilization, and sensitivity across different genera; show the relative ranking of total unique metabolites of a genus in aforementioned categories; give a general overview of the BacDive database stats, and much more. 

For more information about the data available from BacDive, visit the [BacDive Overview](BacDive_Overview.md).

# How to use the streamlit :)
The webapp is split into six tabs: General Overview, Circos, Trends, Cards, By the Numbers, and Comparison.
They each produce different types of visualizations, tables, and/or graphs.

## General Overview
The general overview allows you to do selections for negatively/positively tested, strain/isolate, production/utilization/resistance/sensitivity

You can then select genera by highest absolute values or proportions: that is, genera that have the highest number of species that are ‘resistant’ to the metabolite, or genera that have the highest proportion of species that are. However, we advise the user to take this section with a grain of salt, as it is just meant to provide a quick overview.

## Circos
Users are able to choose a ‘primary’ category of interest — whether the edges be metabolites or genera. They can also choose a ‘secondary’ categor(ies) of interest in order to overlap multiple graphs on a single one. The output will be in .txt files, that you will need to run on circos yourself. You must install it for this to work! Instructions available in [Circos Instructions](Circos_Instructions.md).

## Trends
Two possible visualizations on this page:
1) A parallel diagram. Another way of looking at the number of species in each category (prod/res/sen/util), across genera (that you can choose). 

2) A sankey diagram: select a genus, and you'll see which metabolites it has been tested for, in which category, and whether it was tested positively or negatively.
    
## Cards
Default load 'playing cards' for each of the genera, where you can see how many metabolites the isolate/strain species in the genus have been positively tested for production/utilization/resistance/sensitivity.

## By the Numbers
There are three main features on this page. First is the summary statistics section at the top, which details how many metabolites, genera, isolates, and strains there are in your dataset. Following that, there is a “Metabolite Counts Bar Graph” section, where you can plot the top 15 genera based on neg/pos, isolate/strain, and prod/util/res/sen. This gives an overview of which genera have the most data in that category. See the figure below for an example. The third main feature is the “Homogeneous Metabolite Summary by Genus” section. For a selected metabolite category, this section lists genera where, for at least one metabolite, all species (with a minimum of 5 species) tested uniformly positive or uniformly negative. It also shows which metabolite(s) met that criteria. If a genus appears with both positive and negative results (i.e. mixed), it is omitted. 

## Comparison
This is the section where interesting synergies between genera can be found. 
The user can select up to 10 genera to compare. Upon selecting to include/exclude strains, and the top genera, a table is generated with each genus and its species count as the row name, with columns for which metabolites each genera tested positive for in resistance, sensitivity, production, and utilization. 

More important is the generated shared and synergy summary, which shows the prod/util synergy between the genera, i.e., do any of the genera have species that produce metabolites that the others utilize? Do they have shared utilization, which could indicate that certain metabolites are promoting their growth? Is there any metabolite they are all resistant to?
 
(Note that it is very slow to run, currently. Working on it!)


# If you want to make local changes
make sure requirements are installed :)))
to run, do streamlit run app.py

# Notes to devs:
Notes:
Species Analysis is a little slow to load. Still figuring out why?
