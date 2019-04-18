# burden

## Preparing Data

The burden analysis requires two input files:

  1. **Metadata file** (TSV format)

      Describes the sample in each well of a microplate

  2. **Platereader file** (TSV format)

      Containing the signals (OD600, GFP, et.c) recorded in each well over a time course of grouth

### Metadata file format

This is a tab-separated-values text file (TSV) describing the layout of wells and sample that you must create in a text editor or Excel for each plate that you want to analyze.

Column Descriptions:
* **well** - column (letter) and row (number) of the well on the microplate
* **include** - whether to include this well in the analysis. 0/F/false = no, 1/T/true= yes
* **strain** - the strain/construct that was placed in that well. You _cannot_ use a double underscore "__" in strain names!
** Use **blank** in the **strain** column for wells that contain only medium and no cells.
* **isolate** - number or letter designating different versions of the same strain. For example, different colonies picked after a transformation that might have genetic differences in the plasmids/chromosomes. Different values here represent different _biological replicates_ for a given construct. You do not need to fill in this column for **blanks**.
 **description** - _optional_ column with additional information

Example excerpt from a file showing different types of data and formats:
```text
well  include strain  isolate  description
B1  1  blank    just LB
B2  1  blank    just LB
B3  1  blank    just LB
B4  1  BM1  1 burden monitor strain with empty plasmid, biol repl 1, tech repl 1
B5  1  BM1  1 burden monitor strain with empty plasmid, biol repl 1, tech repl 2
B6  1  BM1  1 burden monitor strain with empty plasmid, biol repl 1, tech repl 3
B7  T  BM1  2 burden monitor with empty plasmid, biol repl 2, tech repl 1
B8  True  BM1  3 burden monitor with empty plasmid, biol repl 2, tech repl 2
B9  True  BM1  3 burden monitor with empty plasmid, biol repl 2, tech repl 3
B10 0  blank  contaminated
B11 0  blank  forgot to inoculate
B12 0  blank  outlier
```
Each of the wells that has the same **strain** and **isolate** values defines a _technical replicate_.

### Platereader file format

This is a tab-separated-values text file (TSV) describing the layout of wells and sample that you must create in a text editor or Excel for each plate that you want to analyze.

IMPORTANT: You may need to run the command `dos2unix` to fix line endings on this file before using it!

Example excerpt from a file showing different types of data and formats:
```text
  A1  A2  A3   ...
0s  0.085  0.0852  0.0867  0.08
249s  0.0849  0.0851  0.0867
496s  0.0849  0.0852  0.0865
744s  0.0849  0.0851  0.0864
992s  0.0849  0.0852  0.0864
1240s  0.0848  0.0852  0.0864
...
52532s  0.0858  0.0863  0.0864
0s  214  210  211
248s  212  211  212
496s  212  211  211
743s  211  209  212
991s  210  209  212
1239s  212  210  211
...
Date of measurement: 2019-02-15/Time of measurement: 18:37:05
...other description lines...
```
Format:
* **1st column** - contains the times in seconds of each measurement. The "s" trailing each number is *optional*."
* **remaining columns** - measurements for a specific well

It is expected that there are 2 or 3 types of measurements recorded in the file in blocks of rows in this order:
* **OD600 absorbance values**
* **GFP fluorescence values**
* **user defined values** _Not fully implemented_

### Platereader file format

## Analyzing Data

The `burden.r` script can be passed arguments to control the analysis of your input files. The most basic command is:

```bash
burden.R --input igem001
```

Try running this on the data found in `examples/igem001`.

## Summary Graphs

This creates a summary graph of growth rate versus GFP production for each strain tested. It can accept multiple input file
prefixes to enable comparing data between multiple `burden.R` runs.

```bash
summary-graph.R --input igem001
```





