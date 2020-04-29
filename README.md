# Burden Assay Protocol 


1. Transform 1 µl of any plasmid of interest into a 50 µl aliquot of chemically competent burden monitor strain cells (CalCl2 chemical competent cells; burden monitor cells made as described in methods section). Incubate transformants overnight at 37C on LB agar plates supplemented with appropriate antibiotics: (1) 50 mg/mL Kanamycin (burden monitor strain is kanamycin-resistant), and (2) the antibiotic(s) in the plasmid backbone. 

3.	Pick an antibiotic-resistant colony from the plate and inoculate in a culture tube with 5 mL LB and appropriate antibiotics. Incubate overnight at 37C (~250rpm). If performing a burden assay the next day, prepare overnight cultures of JEB1204-JEB1208 (BFP control strains, inoculate in 5mL LB supplemented with 50 mg/mL Kanamycin + 20 mg/mL Chloramphenicol).  

    - Important note: cells should be cultured as little as possible (~18 hrs) after picking the colony to avoid mutations, especially with burdensome plasmids. 


4.	Next day: Freeze down 15% glycerol stocks of cultures at -80C and perform burden assay. 

5.	Burden assay (adapted from Ceroni 2015): 
    - **Required materials:**
        - Required samples: overnight cultures of JEB1204 – JEB1208, and any transformed plasmid cultures. 
        - Plate reader (we use Tecan Infinite Pro M200 Plate Reader)  
        - 96 well plate: Black with clear optical bottom and clear lid (we use Nunc™ MicroWell™ 96-Well Optical-Bottom Plates (Thermo Scientific Catalog Number 265301)). 
        - Preheated LB media

### Steps for burden assay:  

 1.	Turn on the plate reader and let it heat to 37C. 
 2.	Create a metadata file (named expXXX.metadata.csv) which will map out the samples in the plate. The file should follow the formatting shown in the image below.  Sample metadata files are in the `examples/input` folder. 
    - XXX denotes unique experiment number/burden assay run  

You will pipette 5µl of each culture in triplicate (for example: one culture is placed in G1, G5, G9, another culture is placed in G2, G6, G10, etc.), as displayed in the diagram below:  
![sample microplate setup](https://github.com/barricklab/burden/blob/master/images/sampleplate.png)   

3. First, prepare the Standard and Blank wells: 
   - Following the diagram shown above, pipette 5µl of each BFP standard (JEB1204-JEB1208) in triplicate (for example: one culture is placed in G1, G5, G9, another culture is placed in G2, G6, G10, etc.).
   - For Blanks (indicated by blue wells): just add 5µl LB  

4. Load 5 µl of each “unique” transformatiion culture in triplicate.
   - This allows up to 23 unique strains to be measured per plate.
   - If you have less than 23 unique strains, load the rest of the wells with 5µl LB, and treat as blanks.  

5. Finally, pipette 195µl preheated LB to all wells, using a multichannel pipette. Pipette contents up and down to ensure they are mixed properly.  

6. Ensure the plate reader has reached 37C, then load the plate in the plate reader.  

7.	Run a burden assay program set to the following parameters: 
    - Record optical density at 600 nm and GFP fluorescence (excitation: 485 nm; emission 528 nm) levels every 10 minutes with 7 minutes of orbital shaking during each cycle. 
    - We use the Magellan software program.
8.	Run each burden assay for a minimum of 6 hours. 
    - If you choose to run for 6 hours (or any amount of time that ends before the program automatically ends), you can “break” the program. Data will not get lost, this action will simply end the run. 

9.	After the run (either after 6 hours or when Magellan/program automatically ends):  

    1.	 An excel sheet should have automatically opened, and it contains the measurements from the completed burden assay. Save this file using the following formatting: “expXXX.measurements.csv” 
    2.	With the Magellan window still open, click “export as ascii” from the file menu, and name the file “expXXX.measurements.asc.” This file will serve as a backup measurements file. 

10.	Create a folder that combines (1) the metadata file from step 2, (2), the measurements file from step 9, and (3), the burden assay program that was used for this particular run (to record the plate reader settings that were used).   

    - Name this folder “expXXX.” (Example: exp057)
      - Ensure that both the metadata and measurements files are saved with the .csv (or .tsv) file extension. 
      - Examples of both file types are found in the `examples/input` folder. 

11.	Next step: Data analysis (computational protocol). 

# Burden analysis for a single data set 

**Required files:** An experiment folder (expXXX) containing the metadata and measurements files from the burden assay.

   - They need be saved with the .csv or .tsv file extension.  
   - (Note that the “XXX” should be a unique burden assay experiment number (such as exp001, exp002, etc). 

1. To analyze a single data set (from one burden assay), download the `burden_fit.R` and  `burden_summary.R` scripts, found in the `scripts` folder. 

   - See the script instructions (https://github.com/barricklab/burden/tree/master/scripts) for a full description of the input file formats and directions for running the scripts at the command line. 

2. Open the `burden_fit.R` file on RStudio. 

      - Download all the required packages: tidyverse, gridExtra, cowplot, optparse.  
 
3.	Set the working directory to the expXXX folder.

4.	Define the variable `input.prefix` in the console, then run the full script. 

      - Example: `input.prefix = “exp057”`

5.	If no errors occur, the expXXX folder should now have many new output files. An example output folder can be found in `examples/exp057`. 

      - *Example output files for exp057 (Anderson series of promoters with RFP)*: 
      
         - `exp057.settings.csv`, a dataframe that prints out all of the options that were used in the fitting.
         - `exp057-plots` folder, which contains three types of files for each strain that was measured: 
            - `samplename.pdf`: A pdf with three plots, one for GFP, one for OD, and one for growth, all as a function of time. These plots were generated using the raw measurement values from the plate reader. 
![sample rates plots output](https://github.com/barricklab/burden/blob/master/images/sample.png)  

            - `samplename.rates.pdf`: A pdf with three plots, each graphing either the calculated ‘specific.growth.rate’, ‘GFP.rate’, or ‘other.rate’ values as a function of time, for each replicate corresponding to the sample:   
![sample rates plots output](https://github.com/barricklab/burden/blob/master/images/sample_rates.png)


         -  `exp057.tidy.metadata.csv` and `exp057.tidy.measurements.csv`- these are the products of the tidyr function, which cleaned up the metadata and measurement files that were initially imported into the script. 
         - `exp057.rates.all.csv`, a file containing rate calculations (growth, GFP, and 'other') for all wells in the assay. The sample file can be found in the `examples/exp057`.
         - `exp057.rates.summary.csv`, which reports the mean rate calculations (of all replicates) for each stain. In this file, each strain has a corresponding mean rate value (GFP, growth, 'other'), with the values of the standard deviations and upper/lower bounds of the 95% confidence intervals of the means. The sample file can be found in the `examples/exp057`.
         **This file will be used to generate the summary graphs for the burden assay (using the `burden_summary.R` script).**  
         
         
         
8. Open the `burden_summary.R` script on RStudio. Ensure the working directory is still set to the correct expXXX folder.

   - Required packages: tidyverse, plotly, gridExtra, cowplot, optparse, xtable. 

9. Define the variable `input.file.string` in the console, then run the script. 
   
   - Example: `input.file.string = “exp057.rates.summary.csv”`.

10. If no errors occur, the expXXX file should now have four new files. The script generates two types of summary graphs: a **growth rate plot**, and a **burden vs growth rate plot**. Both files come in two forms: a pdf and an interactive plotly (html) file. 
   
    - *Example output files for exp057 (Anderson series of promoters with RFP)*: 
   
      - `exp057.rates.summary.growth_rates.pdf` and `exp057.rates.summary.growth_rates.html` - this plots the growth rate of each strain in the assay, with error bars representing the 95% confidence interval of the mean growth rate (mean of all replicates for each strain): 
      
      ![](https://github.com/barricklab/burden/blob/master/images/sample_summary_growth_rates.png)
      
      
      
      - `exp057.rates.summary.burden_vs_growth_rates.pdf` and `exp057.rates.summary.burden_vs_growth_rates.html` - this linear regression plot graphs the mean GFP rate of each strain as a function of mean growth rate. Each point represents a different strain, and has error bars representing the 95% confidence interval for both variables (mean growth and GFP rates):  
      
      ![](https://github.com/barricklab/burden/blob/master/images/sample_burdenvsgrowth.png)
 








