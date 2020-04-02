# covid19
In a parallel git directory run: 
git clone https://github.com/CSSEGISandData/COVID-19.git

Run CoronaVirusPlotter.ipynb in Jupyter Notebook
When you run the cells, it automatically pulls  the lates COVID-19 build then processes it.

If you want the Italian Data:
add a folder called italy and run git clone inside there for the italian data: https://github.com/pcm-dpc/COVID-19
Run: CoronaVirusPlotter_Italy.ipynb

RT-PCR assay analysis:
I also added a new notebook: CoronaVirus_ProbeMap-PCRThermoPredictor.ipynb
This one has all the available PCR primers used throughout the world. It runs these against emerging sequences from GSAID https://www.epicov.org/epi3/frontend#226ddd

It's designed to tell you if a given country's assay can detect all the novel sequences that are emerging in the population. See https://nextstrain.org/ncov to see how these are spreading across the globe.


The code is a little rough, sorry about that.

