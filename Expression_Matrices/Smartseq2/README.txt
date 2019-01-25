Smart-seq2 data of 1787 P. berghei cells that passed quality control.

SS2_pheno.csv contains the pheno or "colData" for each cell:

key column information:

sample_id: same as rownames and colnames in expression data

ShortenedLifeStage: The experimental groups of cells including
bbSpz (bite back/ injected sporozoite)     
EEF (Liver stage)
Merozoite    
oocyst (day 4 oocysts)     
ook (Bolus ookinetes at 18 and 24 hrs)   
ookoo (48 hr ookinetes and oocysts from midgut tissue)   
Ring (Ring culture)
sgSpz (salivary gland Sporozoites)  
Shz (Late stage culture, includes male, female, trophozoites and schizonts)

ShortenedLifeStage2: same as ShortenedLifeStage but with blood stage classifications from SC3 so Shz is split into Male, Female, Trophozoite, Schizont

ShortenedLifeStage3: same as ShortendLifeStage2 but with ookoo split into Ookinete and Oocyst

ShortendLifeStage4: same as ShortendLifeStage3 but with bbSpz and sgSpz combined to just Spz

Colors in fig1 based on ShortendLifeStage2

________________

SS2_counts.csv contains raw counts

SS2_tmmlogcounts.csv contained tmm normalised log counts


