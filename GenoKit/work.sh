GenoKit uORF -d gff.db -g ath.fa -l 6 -r all -f csv -o zhusitao_uORF.csv -s GFF
GenoKit uORF -d gff.db -g ath.fa -l 6 -r all -f csv -o zhusitao_uORF.csv -s GFF -i AT1G01020.6 -m
GenoKit dORF -d gff.db -g ath.fa -l 6 -r all -f csv -o zhusitao_dORF.csv -s GFF -p 8 

GenoKit uorf -d gff.db -g ath.fa -l 6 -r all -f csv -p 8 -o test/zhusitao_uORF.csv -s GFF


time GenoKit uorf -d Araport11_GTF_genes_transposons.Mar202021.gtf -g ath_chr.fa -l 6 -r all -f csv -p 8 -o test/zhusitao_uORF.csv -s GFF

GenoKit uorf -d Araport11_GTF_genes_transposons.Mar202021.gtf -g ath_chr.fa -l 6 -r all -f csv -o zhusitao_uORF.csv -s GFF -i AT1G01020.6 -m

time GenoKit cds -d Araport11_GTF_genes_transposons.Mar202021.gtf -g ath_chr.fa -r all -f csv -p 1 -o test/zhusitao_cds.csv
time GenoKit uorf -d Araport11_GTF_genes_transposons.Mar202021.gtf -g ath_chr.fa -l 1 -r all -f csv -o test/zhusitao_uORF4.csv -s GFF
GenoKit gene -d Araport11_GFF3_genes_transposons.Mar92021.gff -g ath_chr.fa -o test/zhusitao_gene.csv -f csv -s GFF

# 1
time GenoKit create -s gtf -g Araport11_GTF_genes_transposons.Mar202021.gtf -o test/ -p ath 
time GenoKit create -s gff -g Araport11_GFF3_genes_transposons.Mar92021.correct.gff -o test/ -p ath

GenoKit create -s gff -g archetypal.gff -o ./ -p recorrect
GenoKit uorf -d recorrect.gff -g TAIR10.fa -o uORF.csv -f csv -s gff -p 32 -r all -l 6

# cds 
time GenoKit cds -d test/ath.gtf -g ath_chr.fa -r all -f csv -p 1 -o test/zhusitao_cds3.csv
time GenoKit cds -d test/ath.gtf -g ath_chr.fa -r all -f fasta -p 1 -o test/zhusitao_cds3.fa
time GenoKit cds -d test/ath.gtf -g ath_chr.fa -r all -f gff -p 1 -o test/zhusitao_cds3.gff

# transcript
time GenoKit transcript -d test/ath.gtf -g ath_chr.fa -r all -f csv -p 1 -o test/zhusitao_transcript.csv


# promoter 
time GenoKit promoter -d test/ath.gtf -f csv -g ath_chr.fa -l 10 -o test/zhusitao_promoter.csv -u 0 -i AT1G01010 -v 

# terminator
time GenoKit terminator -d test/ath.gtf -f csv -g ath_chr.fa -l 10 -o test/zhusitao_terminator.csv -u 0 -i AT1G01010 -v 

# exon 
time GenoKit exon -d test/ath.gtf -f fasta -g ath_chr.fa -o test/zhusitao_exon.fa -s gff 

# intron 
time GenoKit intron -d test/ath.gtf -f fasta -g ath_chr.fa -o test/zhusitao_intron.fa -s gff

# uorf 
time GenoKit uorf -d test/ath.gtf -g ath_chr.fa -l 1 -r all -f csv -o test/zhusitao_uORF.csv -s gff

# dorf
time GenoKit dorf -d test/ath.gtf -g ath_chr.fa -l 1 -r all -f csv -o test/zhusitao_dorf.csv -s gff
