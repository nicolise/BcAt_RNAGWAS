0.0 Vivian_Bc has Wei Zhang original data from S://
0.1 BotcyninBotcynalideBotrydial.xlsx contains tabs copied from result.lsm.csv (BOT transcripts as phenotypes)
	and XXX (annotations for BOT genes)
	BcBOT_result.lsm.csv 
1.0 bigRRGWAS/01_NameIsolates.R 
	takes IsolateKey_Vivian.csv and BcBOT_result.lsm.csv
	writes BOTphenotypes.csv
2.0 bigRRGWAS/02_MatchGenos.R
	on Linux:
	takes 01_MatchGenos/BOTphenotypes.csv and 01_MatchGenos/hp_bin_trueMAF20_20NA.csv
	hp_bin_trueMAF20_20NA.csv comes from BcSolGWAS/data/GWAS_files/02_csvPrep
	writes 02_bigRR/BOT_phenos.csv and 02_bigRR/BOT_genos.csv
