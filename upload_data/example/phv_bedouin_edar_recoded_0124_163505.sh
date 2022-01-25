#!/bin/bash

echo '# Interpolate physical positions'
echo '---------------------------------------------------------'

/home/acolomer/.conda/envs/phv/bin/tabix -f -p vcf /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/bedouin_edar_recoded.vcf.gz

/home/acolomer/.conda/envs/phv/bin/bcftools query -f '%POS\t%ID\n' /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/bedouin_edar_recoded.vcf.gz > /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/rsid.txt

/home/acolomer/.conda/envs/phv/bin/bcftools query -f '%CHROM\t%ID\t%POS\t%POS\n' /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/bedouin_edar_recoded.vcf.gz > /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/pos_to_interpolate.txt

/usr/bin/Rscript /home/pophumanvar/phv_pipeline/interpolate.R /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/pos_to_interpolate.txt /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/map.txt

if ! [ -s *.vcf.gz.tbi ] | [ -s rsid.txt ] | [ -s pos_to_interpolate.txt ] | [ -s map.txt ];then
	echo '########################################################'
	echo '(!) ERROR! Check your vcf conditions, somthing is wrong'
	echo '########################################################'
fi


echo '# Computing iHS'
echo '---------------------------------------------------------'

/home/pophumanvar/phv_pipeline/selscan-linux-1.3.0/selscan --ihs --vcf /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/bedouin_edar_recoded.vcf.gz --out /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/out --map /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/map.txt --cutoff 0.5 --maf 0.05 --max-gap 200000 --gap-scale 20000 --max-extend 1000000

/home/pophumanvar/phv_pipeline/selscan-linux-1.3.0/norm --ihs --files /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/out.ihs*.out

if ! [ -s out.ihs.out.100bins.norm ];then
	echo '##############################################'
	echo '(!) ERROR! Check iHS warnings in the logfile'
	echo '##############################################'
fi


echo '# Computing nSL'
echo '---------------------------------------------------------'

/home/pophumanvar/phv_pipeline/selscan-linux-1.3.0/selscan --nsl --vcf /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/bedouin_edar_recoded.vcf.gz --out /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/out --maf 0.05 --max-gap 200000 --gap-scale 20000 --max-extend-nsl 100

/home/pophumanvar/phv_pipeline/selscan-linux-1.3.0/norm --nsl --files /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/out.nsl*.out

if ! [ -s out.nsl.out.100bins.norm ];then
	echo '##############################################'
	echo '(!) ERROR! Check nSL warnings in the logfile'
	echo '##############################################'
fi


echo '# Computing iSAFE'
echo '---------------------------------------------------------'

/home/pophumanvar/phv_pipeline/iSAFE/src/isafe.py --input /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/bedouin_edar_recoded.vcf.gz --output /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/out --format vcf --region 2:10950615-12949469 --AA /home/pophumanvar/phv_pipeline/ancestral_alleles/homo_sapiens_ancestor_2.fa --IgnoreGaps --window 300 --step 150 --topk 1 --MaxRank 15 --MaxFreq 0.95

if ! [ -s out.iSAFE.out ];then
	echo '##############################################'
	echo '(!) ERROR! Check iSAFE warnings in the logfile'
	echo '##############################################'
fi


echo '# Preparing email'
echo '---------------------------------------------------------'

if [ -s out.ihs.out.100bins.norm ] && [ -s out.nsl.out.100bins.norm ] && [ -s out.iSAFE.out ];then
	echo '# Merging PHV info'
	/usr/bin/Rscript /home/pophumanvar/phv_pipeline/info.R 2 10950615 12949469 /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/

	if [ -s info.txt.gz ];then
		echo '# Zipping data'
		cp /home/pophumanvar/phv_pipeline/report.rmd /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/Report_chr2-10950615-12949469.rmd
		mkdir /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/data
		cd /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/
		mv $(ls --ignore=log.txt --ignore=*.gz* --ignore=*.rmd* --ignore=*.sh* --ignore=data) /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/data
		cp /home/pophumanvar/phv_pipeline/README ./ 
		echo '# zip information'
		/bin/zip -r Report_chr2-10950615-12949469.zip *
		echo '# Sending results'
		/usr/bin/php /home/pophumanvar/phv_pipeline/mail/phv_mail.php jesus.murga@uab.cat /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/Report_chr2-10950615-12949469.zip
	else
		echo '# Somthing went wrong, sending logfile email'
		/usr/bin/php /home/pophumanvar/phv_pipeline/mail/phv_mail.php jesus.murga@uab.cat /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/log.txt
	fi
else
	echo '# Somthing went wrong, sending logfile email'
	/usr/bin/php /home/pophumanvar/phv_pipeline/mail/phv_mail.php jesus.murga@uab.cat /home/pophumanvar/phv_pipeline/bedouin_edar_recoded_0124_163505/log.txt
fi
