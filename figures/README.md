# Manuscript Figures

## CLK1 Exon 4 splicing visualization
We utilized [ggsashimi tool](https://github.com/guigolab/ggsashimi) to visualize CLK1 exon 4 (see tool repo for detailed instructions). We chose two representative samples BS_30VC3R2Q (high Ex4 inclusion) and BS_FYP2B298 (low Ex4 inclusion).

### Input
```gencode.v39.primary_assembly.annotation.CLK1-select.gtf``` <br>
```input_bams.tsv```

### Generate sashimi plot
1. Clone ggsashimi and pull the official docker image.
```
git clone git@github.com:guigolab/ggsashimi.git
cd ggsashimi
docker pull guigolab/ggsashimi
```

2. Configure the input files as provided in ```ggsashimi/examples```

3. Generate the sashimi plot
```
docker run -w $PWD -v $PWD:$PWD guigolab/ggsashimi -b examples/input_bams.tsv -c chr2:200853000-200864700 -M 15 -C 3 -O 3 --alpha 1 --base-size=10 -R 300 --height=2 --width=8 --ann-height 1 -P examples/palette.txt -o CLK1-ex4-sashimi --shrink -g gencode.v39.primary_assembly.annotation.CLK1-select.gtf
```
