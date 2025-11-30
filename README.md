# Introduction
CNNKSCEC:A Deep Learning–Based Framework for Chromatin Loop Prediction with Multi-Source Feature Integration
# Installation
CNNKSCEC requires Python3 and several scientific packages to run.
```
git clone https://github.com/zhengbingzi/CNNKSCEC.git
conda create -n CNNKSCEC python=3.7.12
conda activate CNNKSCEC
```
To set up a conda virtual environment, the following installation packages are required:
scikit-learn=1.0.2   
tensorflow=2.11.0  
hic-straw=1.3.1
joblib=1.3.1  
numpy  
scipy   
pandas   
h5py   
cooler
pysam  
pyBigWig
tqdm
juicer_tools
# Usage
## Data preparation
Data preparation mainly involves:downloading the.hic file and extracting the hic contact matrix from the.hic file, and generating the hic submatrix;downloading the bigwig file based on DNase-seq and generating the corresponding dnase submatrix.
### Extracting Hi-C contact matrix from.hic file
The process obtains the hic contact matrix for each chromosome from the.hic file.  
This process extracts the hic contact matrix of each chromosome at a specific resolution from the.hic file. It will output frequency_matrix file.
Modify the path to the input and output files in the GetEachchr_Cells_KRobserved.sh file and pass in specific parameters: The.jar file is the path where the juicer tools resides.and run:
```
bash GetEachchr_Cells_KRobserved.sh
```
### Generating sub-matrix from Hi-C contact matrix
The process cuts the Hi-C contact matrix of each chromosome into multiple submatrices.Modify the paths of the input and output files in the Getsubmatrix_chr_all_sample.sh file, where the input is the frequency output from the previous step_matrix file, run:
```
bash Getsubmatrix_chr_all_sample.sh
```
## Model training
If you want to retrain the model, please follow the training data generation method outlined in our paper to obtain the required training samples.
### Get the norm_factor file, Interaction_frequency file
If you don't have norm_factor files and interaction_frequency file, please run the following command. The generated file will be used in the subsequent data generation process. Modify the input and output file paths and specific parameters in GetKRnorm_excerces_factor.sh file and GetKRobserved_excerces.sh respectively, and run:
```
bash GetKRnorm_excerces_factor.sh
bash GetKRobserved_excerces.sh
```
### Get the training positive sample of Hi-C
Here, the files gm12787.tang.ctcf-chiapet.hg19.bedpe and GM12787.mumbach.H3k27ac-hicpip.hg19.bedpe need to be entered. The files are located in the data folder. And run:
```
python get_trainpositive_point.py  -i1 [ctcf.bedpe] -i2 [h3k27ac.bedpe] -o [PositiveTxt_dir] -p [processdata_dir]
python get_trainpositive_sample.py -d [PositiveTxt_dir] -b [frequence_matrix_dir] -o [Positivenpy_dir]  -r [res]
```
The get_trainpositive_point.py file outputs the coordinates of the positive sample center points, while the get_trainpositive_sample.py file outputs the Hi-C positive sample matrix files.
Note that this experiment is conducted at a 5kb resolution. To perform experiments at different resolutions, please modify the resolution value in the corresponding locations in the ctcf_h3k27ac_merge.py and get_trainpositive_point.py files. When running get_trainpositive_sample.py, pass the resolution parameter -r [res] via the command line.
### Get the training negative sample of Hi-C
First, generate different types of negative sample center point files based on the positive sample center point file by running get_trainnegative_point.py, which will output the negative sample center point file.
```
python get_trainnegative_point.py -p [PositiveTxt_dir] -n [norm_factor_dir] -o [NegativeTxt_dir] 
```
Then, split the negative sample center point file for each chromosome into 100 parts for further screening. The valid center points with an interaction frequency greater than 1 are selected, and the filtered negative samples along with their center point files are stored.
```
python get_trainnegative_sample.py -x [frequence_matrix_dir] -p [NegativeTxt_dir] -o [Negativenpy_dir] -n [negative_name_sort_dir] -i [Interaction_frequency_dir] 
```
Merge the negative sample center points for each chromosome by modifying the merge_negpoint_fin.py file. The input file should be the path to the negative_name_sort file. And run:
```
python merge_negpoint_fin.py
```
Merge the negative samples for each chromosome, and run:
```
python merge_chr_npy.py -f [Negativenpy_dir] -o [Negativenpy_Fin_dir]
```
Here,the input file is the path to the Negativenpy file, and the output will be the path to the final negative sample file.
### Get the training positive sample of DNase
Get the DNase-seq matrix corresponding to the positive sample center point file and construct the DNase-seq positive sample file. Modify the GetDnase_Matrix_Positive.py file, with the input file being the positive sample center point file. And run:
```
python GetDnase_Matrix_Positive.py
```
Merge the spatially aligned Hi-C and DNase-seq matrices to create the final training positive sample file.
```
python merge_Dnase_Hic_pos.py
```
### Get the training negative sample of DNase
Get the DNase-seq matrix corresponding to the negative sample center point file and construct the DNase-seq negative sample file. Modify the GetDnase_Matrix_Negative.py file, with the input file being the negative sample center point file. And run:
```
python GetDnase_Matrix_negative.py
```
The space-aligned Hi-C and DNase-seq negative sample files are merged to form the final negative sample file.
```
python merge_Dnase_Hic_neg.py
```
### Generating training and test datasets
Merge Hi-C and DNase positive and negative sample files to generate training and test samples. (In this experiment, training is performed on chromosomes 1-19, and testing is done on chromosomes 20-22).
```
python get_train_test_set.py
```
### Training
```
python train.py
```
### Test
```
python test.py
```
## Prediction
### Get the predicted samples
Extracting the center point files from the Hi-C submatrix of the entire chromosomes generated in the previous work, obtaining the corresponding DNase submatrix based on this file, and merging the Hi-C and DNase submatrices as prediction samples (predictions were made on chromosomes 20-22 in this study). And run:
```
python Getchr_all_center_point.py
python get_all_dnase.py
python merge_Dnase_HiC_all.py
```
### Prediction of chromatin loops using CNNKSCEC
Run the following code to predict the chromatin loops of the entire genome:
```
python cnnkscec_loops.py
```
### Clusting
The clustering method in 
peakachu https://github.com/tariks/peakachu was used for clustering screening:
```
git clone https://github.com/tariks/peakachu
```
Configure the peakachu environment based on the README.md file of peakachu,and run:
```
bash peakachucluster.sh
```
# Output file format
```
[chrname_column1]            The chromosome name1 of the left anchor of the chromatin loop
[location10]                 The starting position of the left anchor of the chromatin loop
[location11]                 The end position of the left anchor of the chromatin loop
[chrname_column2]            The chromosome name2 of the right anchor of the chromatin loop
[location20]                 The starting position of the right anchor of the chromatin loop
[location21]                 The end position of the right anchor of the chromatin loop
[predictions]                The predictions of chromatin loop
[infy]                       The interaction strength of chromatin loop

```




