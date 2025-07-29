# 1
GenoKit create -s gff -g archetypal.gff -o ./ -p recorrect
GenoKit create -s gtf -g chr22.gtf -o ./ -p human
GenoKit create -s gtf -g chr5.gtf -o ./ -p human

# 2
GenoKit uorf -d recorrect.gff -g TAIR10.fa -o uORF.csv -f csv -s gff -p 32 -r all -l 6
GenoKit uorf -d human.gtf -g chr22.fa -o human_uORF_22.csv -f csv -s gtf -p 32 -r all -l 6
GenoKit uorf -d human.gtf -g chr5.fa -o human_uORF_5.csv -f csv -s gtf -p 32 -r all -l 6

# 3
GenoKit sirna -d recorrect.gff -f csv -g chr1.fa -o zhusitao -s gff

# 4
GenoKit vision -i AT1G01210 -d recorrect.gff -l 6 -f pdf -g chr1.fa -o z.pdf -s gff
GenoKit vision -i AT1G01030 -d recorrect.gff -l 6 -f pdf -g chr1.fa -o AT1G01030.pdf -s gff
GenoKit vision -i AT1G05670 -d recorrect.gff -l 6 -f pdf -g chr1.fa -o AT1G05670.pdf -s gff
GenoKit vision -i AT1G06150 -d recorrect.gff -l 6 -f pdf -g chr1.fa -o AT1G06150.pdf -s gff
GenoKit vision -i AT1G80850 -d recorrect.gff -l 6 -f pdf -g chr1.fa -o AT1G80850.pdf -s gff


# SREBP2
GenoKit vision -i ENSG00000198911 -d human.gtf -l 6 -f pdf -g chr22.fa -o ENSG00000198911.pdf -s gtf
# RIMBP3
GenoKit vision -i ENSG00000275793 -d human.gtf -l 6 -f pdf -g chr22.fa -o ENSG00000275793.pdf -s gtf



# HMGCR ENST00000287936.9
GenoKit vision -i ENSG00000113161 -d human.gtf -l 6 -f pdf -g chr5.fa -o ENSG00000113161.pdf -s gtf

# ENSG00000145545
GenoKit vision -i ENSG00000145545 -d human.gtf -l 6 -f pdf -g chr5.fa -o ENSG00000145545.pdf -s gtf

# ENSG00000164318
GenoKit vision -i ENSG00000164318 -d human.gtf -l 6 -f pdf -g chr5.fa -o ENSG00000164318.pdf -s gtf

# ENSG00000142319
GenoKit vision -i ENSG00000142319 -d human.gtf -l 6 -f pdf -g chr5.fa -o ENSG00000142319.pdf -s gtf

# IL3
GenoKit vision -i ENSG00000164399 -d human.gtf -l 6 -f pdf -g chr5.fa -o ENSG00000164399.pdf -s gtf

# CD14
GenoKit vision -i ENSG00000170458 -d human.gtf -l 6 -f pdf -g chr5.fa -o ENSG00000170458.pdf -s gtf

# IRX1
GenoKit vision -i ENSG00000170549 -d human.gtf -l 6 -f pdf -g chr5.fa -o ENSG00000170549.pdf -s gtf



# BRCA1
GenoKit vision -i ENSG00000012048 -d human.gtf -l 6 -f pdf -g chr17.fa -o ENSG00000012048.pdf -s gtf

# 5
GenoKit dorf -d recorrect.gff -g chr1.fa -o dORF.csv -f csv -s gff -p 32 -r all -l 6
GenoKit dorf -d human.gtf -g chr22.fa -o human_dORF_22.csv -f csv -s gtf -p 32 -r all -l 6
GenoKit dorf -d human.gtf -g chr5.fa -o human_dORF_5.csv -f csv -s gtf -p 32 -r all -l 6


# 6
cat chr19.gtf | awk '$3=="transcript"' | awk -F "\t" '{print $9}' | cut -d ";" -f 1| sort | uniq -c | sort -rnk1 | less
