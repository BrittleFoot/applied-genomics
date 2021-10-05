# Deep Sequencing, Error Control


## How 2 Run

1. Download missing srrs and place approprietly
    - `data/SRR1705851.fastq`
    - `data/control/SRR1705858.fastq`
    - `data/control/SRR1705859.fastq`
    - `data/control/SRR1705860.fastq`

2. Install python dependencies
```bash
conda install -c confa-forge -c bioconda -c defaults --file requirements.txt
```

3. Align experimental data to reference genome

```bash
python flu.py prepare data/SRR1705851.fastq data/control
```

4. Do experiments!

> There are couple additional arguments to customize experiment params

```bash
python flu.py compare results results/control/SRR1705858 results/control/SRR1705859 results/control/SRR1705860
```

# Reuslts
I will post results here from the experiment with params freq=0.001 and depth=50000

```
   aa#   results     SRR1705858   SRR1705859   SRR1705860 
    84    0.17%        0.25%        0.19%        0.23%    
   102    0.94%          -            -            -      
   113    0.17%        0.23%        0.21%        0.21%    
   129    0.22%        0.23%          -           0.2%    
   230    0.17%        0.23%        0.23%        0.23%    
   240     0.2%        0.23%        0.23%        0.25%    
   267    0.23%          -            -            -      
   286    0.18%        0.27%        0.29%        0.25%    
   347    0.18%        0.24%          -          0.21%    
   404    0.22%        0.24%        0.22%        0.24%    
   426    0.18%        0.25%        0.24%        0.25%  
```