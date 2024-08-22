#!bin/bash
## If you are working on the Imperial HPC, make sure you follow the conventions of setting up your job script and temporary directories for the job.
## If not, (like I did here, on Silwood Park's William Harvey), proceed with this running this bash script as a job.

#SBATCH --time=0-24:00:00   # Maximum time limit (30 minutes)
#SBATCH --ntasks=1          # Number of tasks
#SBATCH --cpus-per-task=256   # Number of CPU cores per task
#SBATCH --mem=720G            # Memory per node
#SBATCH --partition=large_336

source activate qiime2-amplicon-2024.5

## 3. Import data into QIIME2. For PacBio sequences, the data format is SingleEndFastqManifestPhred33V2. Create a manifest .tsv file, with the first column being sample-id and the second one being absolute-filepath. When creating a manifest .tsv file, ASK COMMAND LINE TO DO IT FOR YOU, SO IT IS EXACTLY THE PATH AS WRITTEN. I made a mistake by copying it in manually, and that took me 30 mins to detect :)
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path $HOME/Data/PB-Loveline2/PB-Loveline/loveline-manifest.tsv \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path $HOME/Results/PB-Loveline2-results/loveline-demux-single-end.qza

## 4. Provide a summary of the imported data
qiime demux summarize \
--i-data $HOME/Results/PB-Loveline2-results/loveline-demux-single-end.qza \
--o-visualization $HOME/Results/PB-Loveline2-results/loveline-demux-single-end-summ.qzv
  
## 5. Trim, filter, denoise, dereplicate, and filter chimeras using DADA2. For more information, go to https://docs.qiime2.org/2022.2/plugins/available/dada2/denoise-ccs/
qiime dada2 denoise-ccs \
  --i-demultiplexed-seqs $HOME/Results/PB-Loveline2-results/loveline-demux-single-end.qza \
  --p-front AGRGTTYGATYMTGGCTCAG \
  --p-adapter RGYTACCTTGTTACGACTT \
  --p-min-len 1000 \
  --p-max-len 1600 \
  --p-chimera-method 'consensus' \
  --o-table $HOME/Results/PB-Loveline2-results/loveline-dada2-results/loveline-feature-table.qza \
  --o-representative-sequences $HOME/Results/PB-Loveline2-results/loveline-dada2-results/loveline-sequence-table.qza \
  --o-denoising-stats $HOME/Results/PB-Loveline2-results/loveline-dada2-results/loveline-dada2-stats.qza \
  --output-dir $HOME/Results/PB-Loveline2-results/loveline-dada2-results \
  --verbose
  
## 6. Visualise DADA2 output using QIIME2View. # Create a metadata .tsv file. For more information, go to https://docs.qiime2.org/2024.5/tutorials/metadata/. 
## Essentially, do it in Google Sheets - there is a special tool that helps check your metadata file, to ensure QIIME2 can work with it. 
qiime metadata tabulate \
  --m-input-file $HOME/Results/PB-Loveline2-results/loveline-dada2-results/loveline-dada2-stats.qza \
  --o-visualization $HOME/Results/PB-Loveline2-results/loveline-dada2-results/dada2-ccs_stats.qzv

qiime feature-table summarize \
  --i-table $HOME/Results/PB-Loveline2-results/loveline-dada2-results/loveline-feature-table.qza \
  --o-visualization $HOME/Results/PB-Loveline2-results/loveline-dada2-results/loveline-feature-table.qzv \
  --m-sample-metadata-file $HOME/Data/PB-Loveline2/PB-Loveline/loveline-metadata.tsv
  
qiime feature-table tabulate-seqs \
  --i-data $HOME/Results/PB-Loveline2-results/loveline-dada2-results/loveline-sequence-table.qza \
  --o-visualization $HOME/Results/PB-Loveline2-results/loveline-dada2-results/loveline-sequence-table.qzv \
  --output-dir loveline-dada2-results-table
  
## 7. Cluster sequences de novo to create an OTU table using VSEARCH
qiime vsearch cluster-features-de-novo \
  --i-sequences $HOME/Results/PB-Loveline2-results/loveline-dada2-results/loveline-sequence-table.qza \
  --i-table $HOME/Results/PB-Loveline2-results/loveline-dada2-results/loveline-feature-table.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table $HOME/Results/PB-Loveline2-results/loveline-cluster-results/loveline-clust-feature-table.qza \
  --o-clustered-sequences $HOME/Results/PB-Loveline2-results/loveline-cluster-results/loveline-clust-sequence-table.qza \
  --output-dir $HOME/Results/PB-Loveline2-results/loveline-cluster-results \
  --verbose

## 8. Assign taxonomy using VSEARCH, using the Silva full length classifier from https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza
### 8a. Get the classifier from https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza. More classifiers can be found on https://resources.qiime2.org/.
qiime feature-classifier classify-consensus-vsearch \
 --i-query $HOME/Results/PB-Loveline2-results/loveline-cluster-results/loveline-clust-sequence-table.qza \
 --i-reference-reads $HOME/Data/PB-Loveline2/PB-Loveline/silva-138-99-seqs.qza \
 --i-reference-taxonomy $HOME/Data/PB-Loveline2/PB-Loveline/silva-138-99-tax.qza  \
 --p-threads 256 \
 --p-maxrejects 100 --p-maxaccepts 10 --p-perc-identity 0.97 \
 --o-classification $HOME/Results/PB-Loveline2-results/loveline-tax-table.qza

## 9. Convert the clust-feature-table.qza to a presence-absence table. We do this by taking the clust-sequence-table.qza
## we generated beforehand and making a relative frequency table before converting to a presence absence table.
qiime feature-table relative-frequency \
  --i-table $HOME/Results/PB-Loveline2-results/loveline-cluster-results/loveline-clust-feature-table.qza \
  --o-relative-frequency-table $HOME/Results/PB-Loveline2-results/loveline-clustseq-table-relfreq.qza

qiime feature-table presence-absence \
  --i-table $HOME/Results/PB-Loveline2-results/loveline-clustseq-table-relfreq.qza \
  --o-presence-absence-table $HOME/Results/PB-Loveline2-results/loveline-clustseq-table-presabs.qza

qiime taxa barplot \
  --i-table $HOME/Results/PB-Loveline2-results/loveline-clustseq-table-presabs.qza \
  --i-taxonomy $HOME/Results/PB-Loveline2-results/loveline-tax-table.qza \
  --o-visualization $HOME/Results/PB-Loveline2-results/loveline-tax-table.qzv
  
