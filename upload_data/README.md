<h2>Upload tab</h2>
<b>PopHumanVar</b> provides and option to process the users data.<br>
In the scripts seciton you can find some of the scripts used in the pipeline.<br>
In the exemple directory you can find an example of the results the user gets on their email if everything works fine. It incluedes:

- <b>README:</b> explanation of how the data is organized in the directory, and how to proceed to see the results.
- <b>vcf.gz + vcf.gz.tbi</b>: vcf file uploaded (not included here)
- <b>info.txt.gz</b>: compressed file with the informaiton from all ddbb
- <b>phv_filename.sh</b>: pipeline
- <b>log.txt</b>: summary of the pipeline and errors encountered
- <b>Report_region.rmd</b>: simplified copy of PopHumanVar for Rstudio
- <b>data:</b> logs and files generated in intermediate steps 
    - <b>pos_to_interpolate.txt</b>: list of	chr, rsid and physical position (intermediat step for interpolation)
    - <b>rsid.txt</b>: physical position & rsid of each entry in your data
    - <b>map.txt</b>: recombination map after interpolation
    - <b>out.ihs.out</b>: raw results from selscan (ihs)
    - <b>out.ihs.log</b>: report of selscan process (ihs)
    - <b>out.ihs.out.100bins.norm</b>: normalitzation of ihs results above
    - <b>out.nsl.out</b>: raw results from selscan (nsl)
    - <b>out.nsl.log</b>: report of selscan process (nsl)
    - <b>out.nsl.out.100bins.norm</b>: normalization of the nsl results above
    - <b>out.iSAFE.out</b>: results	from iSAFE
    - <b>logfile</b>: report of isafe process
<br>
It was done using the <a href="https://pophumanscan.uab.cat/data/phv/bedouin_edar_recoded.vcf.gz">file</a> provided as an example in the app.<br>
For more information, check the <a href="https://pophumanvar.uab.cat/#tutorial">tutorial</a> from PopHumanVar.
