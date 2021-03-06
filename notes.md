# Notes for Calakmul Landsat Time Series Analysis

## Important references

The paper from DeVries et al. (2016) "Characterizing forest change..." is useful to: a) copy part of the methods (especially pre-processing), and b) see how surface reflectance images are transformed via tasseled cap (Crist's "A Tasseled Cap equivalent transformation for reflectance factor data" (1985) is referenced for the coefficients). BFAST here is used to monitor change, and as such, the study is designed to detect a single disturbance - in other words, they use BFAST in its "online" mode.

The difference between offline (my method) and online modes of change detection is described (briefly) in Zhu's review "Change detection using Landsat time series" (2017).

The Distrubance Index, derived from tasseled cap bands, is described in Healey et al's "Comparison of Tasseled Cap-based Landsat data structures..." (2005). In the paper they argue that tasseled cap-derived data is better to detect forest disturbances than original Landsat data. Accuracy did not var greatly among tasseled cap data structures. Note that to use DI, the scene needs to be standarized - something that might cause problems in cases where gaps are considerable. Another option is to use the Tasseled Cap Angle (see Powell et al 2010 "Quantification of live aboveground forest biomass dynamics...").

Refer to PengSong et al (2019) "Global land change from 1982 to 2016" to see how Mann-Kendall tests and Sen-Theil slopes are used in the context of forest cover change.


## Reflections for writing

I need to describe local land-use patterns in order to illustrate the relatively simple land conversions taking place in Calakmul, i.e. from forest (including secondary) to agriculture (including pasture) and back. Other conversions can be considered virtually absent in this region. This description is to justify the simplistic approach followed to train the model: to determine prescence or absence of forest cover in a given pixel *visually* (see Cohen's paper on their "TimeSync"). 


## Google Cloud

### R and RStudio Server resources
https://grantmcdermott.com/2017/05/30/rstudio-server-compute-engine/ (setup through gcloud, then launches RStudio Server)

https://cloudyr.github.io/googleComputeEngineR/index.html (everything done from RStudio (desktop) to then launch RStudio Server)

### Add storage
https://cloud.google.com/compute/docs/disks/add-persistent-disk (add, format, and mount instructions)

### Pricing
[Virtual machines](https://cloud.google.com/compute/vm-instance-pricing)  
[Persistent disks](https://cloud.google.com/compute/disks-image-pricing)  

### SDK documentation
https://cloud.google.com/sdk/docs/

### Bash commands
ls: list files
df -h: show mounted drives
 



### G Cloud Commands 
[Link](https://cloud.google.com/sdk/docs/quickstart-debian-ubuntu)to install G Cloud SDK.  
If these instructions were followed, then just type:  
`gcloud init`  

Otherwise, after creating the instance and installing G Cloud SDK:  
`sudo gcloud init --console-only`  > 1 > Y > copy-paste link in browser > copy verification code and paste in terminal > pick project  

SSH into the instance:  
`sudo gcloud compute ssh ubuntu-01 --zone us-west2-a`  
(pph: sethgorelick)  

Mount additional drive  
`sudo mount -o discard,defaults /dev/sdb /mnt/disks/input_data`  

Transfer files between physical computer and instance:  
Physical -> Instance (single file):  
`sudo gcloud compute scp file/in/physical/xxx.RDS name_of_instance:/destination/path --zone zone` (omit zone if pre-established in init)  
e.g. `sudo gcloud compute scp model_randfor_full.RDS ubuntu-01:/mnt/disks/input_data --zone us-west2-a`  

Physical -> Instance (directory):
`sudo gcloud compute scp --recurse folder_name/ name_of_instance:/destination/path --zone zone` (omit zone if pre-established in init)  
e.g. `sudo gcloud compute scp --recurse 00_cropped_masked_bands/ ubuntu-01:/mnt/disks/input_data --zone us-west2-a`

Instance -> Physical (single file):
`sudo gcloud compute scp name_of_instance:/file/in/instance/xxx.RDS /destination/path/in/physical/ --zone us-west2-a` (omit zone if pre-established in init)  
e.g. `sudo gcloud compute scp ubuntu-01:/mnt/disks/input_data/tbb.RDS /media/cdobler/Neobxbaumia/Research/Project_calakmul_landsat_ts/out_training/ --zone us-west2-a`


### Questions to Seth
RStudio setup  
Disaggregate storage costs from instance/computing costs  







