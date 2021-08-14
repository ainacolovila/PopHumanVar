#!/bin/bash 
 
user='shiny'
pass='***'
db=functional
host=localhost
    
commands="CREATE DATABASE ${db};GRANT USAGE ON *.* TO '${user}'@'${host}' IDENTIFIED BY '${pass}';GRANT ALL privileges ON \`${db}\`.* TO '${user}'@'${host}';FLUSH PRIVILEGES;"

echo "${commands}" | /usr/bin/mysql -u root -p

#######################Statistics
for i in `seq 1 22`;
do
    echo ${i}
    cmd1="CREATE table chr${i}(physicalPos INT,rsid VARCHAR(255),ACB FLOAT NULL,ASW FLOAT NULL,BEB FLOAT NULL,CDX FLOAT NULL,CEU FLOAT NULL,CHB FLOAT NULL,CHS FLOAT NULL,ESN FLOAT NULL,FIN FLOAT NULL,GBR FLOAT NULL,GIH FLOAT NULL,GWD FLOAT NULL,IBS FLOAT NULL,ITU FLOAT NULL,JPT FLOAT NULL,KHV FLOAT NULL,LWK FLOAT NULL,MSL FLOAT NULL,PJL FLOAT NULL,STU FLOAT NULL,TSI FLOAT NULL,YRI FLOAT NULL,primary key(physicalPos));"

    cmd2="LOAD data local infile '/home/selector/shiny-server/ShinyBeta/data/${db}/${db}_chr${i}.txt' into table chr${i};"
    
    mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd1}"

    mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd2}"
done

#######################Geva

for i in `seq 1 21`;
do
    echo ${i}
    sed -i '1d' /home/selector/shiny-server/ShinyBeta/data/atlas_variant_age/atlas_chr${i}.txt
    cmd1="CREATE table chr${i}(physicalPos INT,rsid VARCHAR(255),AlleleRef VARCHAR(255),AlleleAlt VARCHAR(255),AlleleAnc VARCHAR(255),AgeMode_Mut FLOAT NULL,AgeMean_Mut FLOAT NULL,AgeMedian_Mut FLOAT NULL,AgeCI95Lower_Mut FLOAT NULL,AgeCI95Upper_Mut FLOAT NULL,AgeMode_Rec FLOAT NULL,AgeMean_Rec FLOAT NULL,AgeMedian_Rec FLOAT NULL,AgeCI95Lower_Rec FLOAT NULL,AgeCI95Upper_Rec FLOAT NULL,AgeMode_Jnt FLOAT NULL,AgeMean_Jnt FLOAT NULL,AgeMedian_Jnt FLOAT NULL,AgeCI95Lower_Jnt FLOAT NULL,AgeCI95Upper_Jnt FLOAT NULL,QualScore_Mut FLOAT NULL,QualScore_Rec FLOAT NULL,QualScore_Jnt FLOAT NULL,key(physicalPos));"

    cmd2="LOAD data local infile '/home/selector/shiny-server/ShinyBeta/data/atlas_variant_age/atlas_chr${i}.txt' into table chr${i};"
    
    mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd1}"
    
    mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd2}"

done
#######################PHS

for i in `seq 1 22`;
do
    echo ${i}
    sed -i '1d' /home/selector/phs/phs_chr${i}.txt

    cmd1="CREATE table chr${i}(start INT NULL,end INT NULL,geneid VARCHAR(255) NULL,pops VARCHAR(255) NULL,source VARCHAR(255) NULL,key(start,end));"

    cmd2="LOAD data local infile '/home/selector/phs/phs_chr${i}.txt' into table chr${i};"
    
    mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd1}"
    
    mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd2}"

done

#######################FUNCTIONAL INFO

for i in `seq 2 22`;
do
    echo ${i}
    sed -i '1d' /home/selector/functional/functional_phv_${i}.txt
    
    cmd1="CREATE table chr${i}(physicalPos INT,rsid VARCHAR(255) NULL,gene_effect MEDIUMTEXT NULL, effect MEDIUMTEXT NULL, impact MEDIUMTEXT NULL, TFbinding INT NULL, DNasePeak INT NULL, motif INT NULL,DNaseFootprint INT NULL,eQTL INT NULL,matchedTFmotif INT NULL,matchedDNaseFootprint INT NULL,ranking VARCHAR(255) NULL,risk_allele VARCHAR(255) NULL,risk_allele_freq TEXT NULL,trait TEXT NULL,accession TEXT NULL,context VARCHAR(255) NULL,diseaseName MEDIUMTEXT NULL,diseaseId TEXT NULL,diseaseType VARCHAR(255) NULL,source VARCHAR(255) NULL,DSI INT NULL,DPI INT NULL,EI INT NULL,NofPmids INT NULL,ClinicalSignificance VARCHAR(255) NULL,ClinSigSimple VARCHAR(255) NULL,PhenotypeList TEXT NULL,key(physicalPos));"

    cmd2="LOAD data local infile '/home/selector/functional/functional_phv_${i}.txt' into table chr${i};"
    
    mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd1}"
    
    mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd2}"

done

#######################SLIDER


cmd1="CREATE table sliderStaticValues(chr INT,variable VARCHAR(255) NULL,minV INT NULL, maxV INT NULL, key(chr));"

cmd2="LOAD data local infile '/home/selector/sliderStaticValues.tsv' into table sliderStaticValues;"

mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd1}"

mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd2}"


#######################CUTOFF   

cmd1="CREATE table statCutoff(
 test VARCHAR(255) NULL, ACB FLOAT NULL, ASW FLOAT NULL, BEB FLOAT NULL, CDX FLOAT NULL, CEU FLOAT NULL, CHB FLOAT NULL, CHS FLOAT NULL, ESN FLOAT NULL, FIN FLOAT NULL, GBR FLOAT NULL, GIH FLOAT NULL, GWD FLOAT NULL, IBS FLOAT NULL, ITU FLOAT NULL, JPT FLOAT NULL, KHV FLOAT NULL, LWK FLOAT NULL, MSL FLOAT NULL, PJL FLOAT NULL, STU FLOAT NULL, TSI FLOAT NULL, YRI FLOAT NULL,key(test));"

cmd2="LOAD data local infile '/home/selector/stat_cutoff.txt' into table statCutoff;"

mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd1}"

mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd2}"


#######################GWAS

# BEFORE PUSH
# sed -i '1d' gwas_catalog.txt 

# sed 's/NA;//g' gwas_catalog.txt  | sed 's/;NA//g' > tmp && mv tmp gwas_catalog.txt 

cmd1="CREATE table gwas_catalog(chr INT,physicalPos INT,rsid VARCHAR(255) NULL,risk_allele VARCHAR(255) NULL,risk_allele_freq VARCHAR(255) NULL,trait  VARCHAR(255) NULL,accession VARCHAR(255) NULL,context VARCHAR(255) NULL,key(chr,physicalPos));"

cmd2="LOAD data local infile '/home/selector/shiny-server/ShinyBeta/data/gwas_catalog.txt' into table gwas_catalog;"

mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd1}"

mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd2}"



#######################DISGENET


cmd1="CREATE table disgenet(chr INT,physicalPos INT,rsid VARCHAR(255) NULL,diseaseName VARCHAR(255) NULL,diseaseId VARCHAR(255) NULL,diseaseType VARCHAR(255) NULL,source VARCHAR(255) NULL,DSI FLOAT NULL,DPI FLOAT NULL,EI  FLOAT NULL,NofPmids INT NULL,key(chr,physicalPos));"

cmd2="LOAD data local infile '/home/selector/shiny-server/ShinyBeta/data/disgenet.txt' into table disgenet;"

mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd1}"

mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd2}"

#######################CLINVAR

cmd1="CREATE table clinvar(chr INT,physicalPos INT,rsid VARCHAR(255) NULL,ClinicalSignificance VARCHAR(255) NULL,ClinSigSimple VARCHAR(255) NULL,PhenotypeList VARCHAR(255) NULL,key(chr,physicalPos));"

cmd2="LOAD data local infile '/home/selector/shiny-server/ShinyBeta/data/Final_Aina_Collapse_Clinvar.txt' into table clinvar;"

mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd1}"

mysql --user="${user}" --password="${pass}" --database="${db}" --execute="USE ${db};${cmd2}"

