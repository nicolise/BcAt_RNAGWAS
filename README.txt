0.0 Vivian_Bc has Wei Zhang original data from S://
0.1 BotcyninBotcynalideBotrydial.xlsx contains tabs copied from result.lsm.csv (BOT transcripts as phenotypes)
	and XXX (annotations for BOT genes)
	BcBOT_result.lsm.csv 
1.0 bigRRGWAS/01_NameIsolates.R 
	takes IsolateKey_Vivian.csv and BcBOT_result.lsm.csv
	writes BOTphenotypes.csv
1.1 bigRRGWAS/02_PrepGenos.R
	tab files from BcSolGWAS/data/GWAS_files/01_tabfiles/Suzi_033016/
	copied to data/00_tabfiles/Suzi_033016/
	output: 02_MatchGenos/
2.0 bigRRGWAS/03_MatchGenos.R
	on Linux:
	takes 02_MatchGenos/BOTphenotypes.csv and 02_MatchGenos/hp_bin_trueMAF20_20NA.csv
	Key_SNPnames.csv comes from BcSolGWAS/data/GWAS_files/02_csvPrep/
	writes 03_bigRR/BOT_phenos.csv and 03_bigRR/BOT_genos.csv

