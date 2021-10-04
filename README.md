# CircSMG
Visualization of strain specific genomic regions and their relation with metabolic functions 

# Requirements

- Linux (not tested on Windows and MacOS)
- [Phyton](https://www.python.org/) 3.8 or later

## The required packages 

- [Pandas](https://github.com/pandas-dev/pandas) pandas-1.3.3
- [Circos](http://circos.ca/) v 0.69-8 


# Installation 

I recommend to use conda but not necessary

```
conda update conda --all
conda create --name CircSMG python=3.8
conda activate CircSMG

```

You can use pip to install dependecies

```
pip install pandas
```
To install circos please visit the site: [Circos] (http://circos.ca/)


# Demo Run

You can try the initially prepared config file and input files to see what can be done with CircSMG. All you need is 
just to run CircSMG.py file to generate Circos inputs.
```
conda activate CircSMG
git clone https://github.com/recepcanaltinbag/CircSMG.git
cd CircSMG
python3 CircSMG.py
```
Then you need to run Circos to create output figure:

```
cd Circos
circos -conf circos.conf
```
You can see the output file in the Circos folder:


![example_output](/Circos/circos.png)
