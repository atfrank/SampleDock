---
layout: default
---

For installation, please follow the instructions in the github [repo](https://github.com/atfrank/SampleDock).<br>
Sample and Dock webservice can be accessed from [SMALTR](https://smaltr.org/).

## Sample and Dock Design Visualizations

### CDK2 nDesign TMAP
Comparison of different number of designs per cycle for CDK2 designs. Only best designs of each cycle is plotted.
<div>
    <iframe src="./vis_maps/tmap_CDK2_ndesigns.html" style="height:800px;width:800px;" title="TMAP of Varing Design per Cycle for CDK2"></iframe> 
    <p style="text-align:right">
        <a href="https://atfrank.github.io/SampleDock/vis_maps/tmap_CDK2_ndesigns.html">View in fullscreen</a> 
    </p>
</div>

### M<sup>Pro</sup> TMAP 
Best design of each cycle for M<sup>Pro</sup> (20 designs each cycle and 24-hour run)
<div>
    <iframe src="./vis_maps/tmap_mpro.html" style="height:800px;width:800px;" title="TMAP of MPro Designs"></iframe> 
    <p style="text-align:right">
        <a href="https://atfrank.github.io/SampleDock/vis_maps/tmap_mpro.html">View in fullscreen</a> 
    </p>
</div>

### CDK2 20 initial Scaffolds TMAP 
Comparison of designs from 20 different initial scaffolds. All designs from each run are plotted (20 designs each cycle and 24-hour run). The 20 scaffolds used are shown below.
<div>
    <iframe src="./vis_maps/tmap_CDK2_20_Scaffolds.html" style="height:800px;width:800px;" title="TMAP of 20 Different Initial Scaffolds for CDK2"></iframe> 
    <p style="text-align:right">
        <a href="https://atfrank.github.io/SampleDock/vis_maps/tmap_CDK2_20_Scaffolds.html">View in fullscreen</a> 
    </p>
    <p style="text-align:center">The Top 20 Scaffolds in ChEMBL Used as Initial Scaffold</p>
    <img src="./vis_maps/20-scaff-wide.svg" alt="Top 20 Scaffold on ChEMBLt" width="800"> 
</div>

------------------------------------------
### Other TMAPS
#### [CDK2 20 Design TMAP](https://atfrank.github.io/SampleDock/vis_maps/tmap_CDK2_20designs_all.html)
All designs for CDK2 from the trial with 20 designs per cycle and 24-hour run

#### [CDK2 Ensemble TMAP](https://atfrank.github.io/SampleDock/vis_maps/tmap_CDK2_ensemble.html)
Best 5 designs of each cycle for CDK2 using top 5 ensembling appproach (20 designs generated from previous top 5 designs for each cycle, and 24-hour run)

#### [CDK2 Latent Space UMAP](https://atfrank.github.io/SampleDock/vis_maps/umap_CDK2_ndesigns.html) 
Comparison of different number of designs per cycle for CDK2 designs. 3D UMAP generated using the latent space coordinates (Tree Space, 28 dimensions). All designs plotted.<br>

Visualization Powered by [Faerun](https://github.com/reymond-group/faerun-python)