# burden

## Preparing Data

The burden analysis requires two input files:

  1. **Metadata file** (CSV or TSV format)

      Describes the sample in each well of a microplate

  2. **Platereader file** (CSV or TSV format)

      Containing the signals (OD600, GFP, et.c) recorded in each well over a time course of grouth

### Metadata file format

This is a comma-separated values (CSV) or tab-separated values text file (TSV) describing the layout of wells and sample that you must create in a text editor or Excel for each plate that you want to analyze.

Column Descriptions:
* **well** - column (letter) and row (number) of the well on the microplate
* **strain** - the strain/construct that was placed in that well. You _cannot_ use a double underscore "__" in strain names!
   * Use **blank** or *b* in the **strain** column for wells that contain only medium and no cells.
* **include** - _optional_ column determining whether to include this well in the analysis. Set it to false to remove wells with technical problems. (0/F/false/blank = no, 1/T/true= yes)
* **isolate** - _optional_ column with a number, letter, or name designating different versions of the same strain. For example, different colonies picked after a transformation that might have genetic differences in the plasmids/chromosomes could be labeled A, B, C,. Different values here represent different _biological replicates_ for a given construct. This value will be ignored for **blanks**.
* **description** - _optional_ column with additional information

Example excerpt from a file showing different types of data and formats:
```text
well  strain  control include isolate description
B1    blank   F       T               just LB
B2    blank   F       T               just LB
B3    blank   F       T               just LB
B4    Test1   F       T       A       test strain 1, transformant A, tech repl 1
B5    Test1   F       T       A       test strain 1, transformant A, tech repl 2
B6    Test1   F       T       A       test strain 1, transformant A, tech repl 3
B7    Test2   F       T       B       test strain 1, biol repl 2, tech repl 1
B8    Test2   F       T       B       test strain 1, biol repl 2, tech repl 2
B9    Test2   F       T       B       test strain 1, biol repl 2, tech repl 3
B10   blank   F       F               contaminated
B11   blank   F       F               forgot to inoculate
B12   blank   F       F               outlier
C1    Ctrl1   T       T               control strain 1, tech repl 1
C2    Ctrl1   T       T               control strain 1, tech repl 2
C3    Ctrl1   T       T               control strain 1, tech repl 3
C4    Ctrl2   T       T               control strain 2, tech repl 1
C5    Ctrl2   T       T               control strain 2, tech repl 2
C6    Ctrl2   T       T               control strain 2, tech repl 3
```
Each of the wells that has the same **strain** and **isolate** are considered together as _technical replicates_.

### Platereader file format

This is a comma-separated values (CSV) or tab-separated values text file (TSV) describing the layout of wells and sample that you must create for each plate that you want to analyze.

IMPORTANT: You may need to run the command `dos2unix` to fix line endings on this file before using it!

Example excerpt from a file showing different types of data and formats:
```text
time   A1      A2      A3      ...
0s     0.0850  0.0852  0.0867  ...
249s   0.0849  0.0851  0.0867  ...
496s   0.0849  0.0852  0.0865  ...
744s   0.0849  0.0851  0.0864  ...
992s   0.0849  0.0852  0.0864  ...
1240s  0.0848  0.0852  0.0864  ...
...
52532s 0.0858  0.0863  0.0864  ...
0s        214     210     211  ...
248s      212     211     212  ...
496s      212     211     211  ...
743s      211     209     212  ...
991s      210     209     212  ...
1239s     212     210     211  ...
...
Date of measurement: 2019-02-15/Time of measurement: 18:37:05
...other description lines...
```
Format:
* **optional header row** – if missing, assigns the first column to time and then 96 wells in the order: A1, A2, A3...
* **1st column** - contains the times *in seconds* of each measurement. The "s" trailing each number is *optional*."
* **remaining columns** - measurements for a specific well*


It is expected that there are 2 or 3 types of measurements recorded in the file in blocks of rows in this order:
* **OD600 absorbance values**
* **GFP fluorescence values**
* **user defined values (like BFP or RFP)**

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
summary-graph.R --input igem001.rates.summary.csv
```

It can accept multiple input file
prefixes to enable comparing data between multiple `burden.R` runs if you separate them with commas (no spaces allowed by the comma)

```bash
summary-graph.R --input "igem001.rates.summary.csv,igem002.rates.summary.csv"
```




