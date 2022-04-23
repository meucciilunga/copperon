/home/katanga/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/tblastn \
-query '/home/katanga/Coding/copperon/blast/queries/CopY_query.faa' \
-db prok_genome_db \
-max_target_seqs 25000 \
-outfmt "6 sseqid sstart send length sframe pident evalue qseqid qstart qend" \
-num_threads 24 \
-db_gencode 11 \
-evalue 1E-12 \
-out '/home/katanga/Coding/copperon/blast/results/CopY_blast_result.txt' &&


/home/katanga/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/tblastn \
-query '/home/katanga/Coding/copperon/blast/queries/CupA_query.faa' \
-db prok_genome_db \
-max_target_seqs 25000 \
-outfmt "6 sseqid sstart send length sframe pident evalue qseqid qstart qend" \
-num_threads 24 \
-db_gencode 11 \
-evalue 1E-12 \
-out '/home/katanga/Coding/copperon/blast/results/CupA_blast_result.txt' &&


/home/katanga/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/tblastn \
-query '/home/katanga/Coding/copperon/blast/queries/CopA_query.faa' \
-db prok_genome_db \
-max_target_seqs 25000 \
-outfmt "6 sseqid sstart send length sframe pident evalue qseqid qstart qend" \
-num_threads 24 \
-db_gencode 11 \
-evalue 1E-12 \
-out '/home/katanga/Coding/copperon/blast/results/CopA_blast_result.txt' &&


/home/katanga/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/tblastn \
-query '/home/katanga/Coding/copperon/blast/queries/CopZ-TcrZ_query.faa' \
-db prok_genome_db \
-max_target_seqs 25000 \
-outfmt "6 sseqid sstart send length sframe pident evalue qseqid qstart qend" \
-num_threads 24 \
-db_gencode 11 \
-evalue 1E-9 \
-out '/home/katanga/Coding/copperon/blast/results/CopZ-TcrZ_blast_result.txt'