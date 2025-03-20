# Characterization of aberrant splicing in pediatric central nervous system tumors reveals CLK1 as a candidate oncogenic dependency
Ammar S. Naqvi, Ryan J. Corbett, Priyanka Seghal, Karina L. Conkrite, Komal S. Rathi, Brian M. Ennis, Katharina E Hayer, Bo Zhang, Miguel A. Brown, Daniel P. Miller, Adam A. Kraya, Joseph M. Dybas, Zhuangzhuang Geng, Christopher Blackden, Shebheel Arif, Antonia Chroni, Aditya Lahiri, Madison L. Hollawell, Phillip B. Storm, Jessica B. Foster, Matuesz Koptyra, Peter J. Madsen, Sharon J. Diskin, Andrei Thomas Tikhonenko, Adam C. Resnick, Jo Lynne Rokita

This project originated at the Center for Data-Driven Discovery in Biomedicine at Children's Hospital of Philadelphia.
For issue and pull request history, please see https://github.com/d3b-center/pbta-splicing.

### Clone repository
```
git clone git@github.com:rokitalab/clk1-splicing.git
```

### Download data
```
bash download_data.sh
```

### Docker set-up
#### Pull docker image
```
docker pull pgc-images.sbgenomics.com/rokita-lab/pbta-splicing:v1.0.0
```
#### docker run
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/clk1-splicing pgc-images.sbgenomics.com/rokita-lab/pbta-splicing:v1.0.0
```
#### docker execute
```
docker exec -ti clk1-splicing bash
```

### Generate paper figures
```
bash scripts/run_analyses.sh
```

### Code Authors
Ammar S. Naqvi (@naqvia), Jo Lynne Rokita (@jharenza), Ryan Corbett (@rjcorb)

### Contact
For questions, please submit an issue.
