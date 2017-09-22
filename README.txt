0.0 Vivian_Bc has Wei Zhang original data from S://
0.1 BotcyninBotcynalideBotrydial.xlsx contains tabs copied from result.lsm.csv (BOT transcripts as phenotypes)
	and XXX (annotations for BOT genes)
	BcBOT_result.lsm.csv 
1.0 bigRRGWAS/01_NameIsolates.R 
	takes IsolateKey_Vivian.csv and BcBOT_result.lsm.csv
	writes BOTphenotypes.csv
2.0 bigRRGWAS/02_MatchGenos.R
	takes BOTphenotypes.csv and
	writes