FILENAME algae '/folders/myfolders/Project1/algae.csv'; /*create a pointer to data file*/
 
DATA algae; /*Assign name algae to data*/
INFILE algae DSD FIRSTOBS = 2;
INPUT season riversize fluidvelocity chem1 chem2 chem3 chem4 chem5 chem6 chem7 chem8 abundance;
RUN;

PROC PRINT DATA=algae (OBS=10); /* Print 10 observations from the original dataset */
TITLE 'Algae dataset';
run;

/*Part a) Finding associaton between riversize and season*/
PROC FREQ DATA=algae;
TABLES riversize*season / CHISQ FISHER; /* contigency table of riversize by season and chisquare test */
TITLE 'Association between riversize and season';
RUN;

/*Part b) Creating a new variable by combining the small and medium size rivers in one category */
DATA algae1; SET algae;
IF riversize= 1 OR riversize= 2 THEN group = 0;
ELSE group = riversize;
RUN;

PROC PRINT DATA=algae1 (OBS=10); /* Print 10 observations from the new dataset */
TITLE 'Updated algae dataset: combining small and medium size rivers in one category';
RUN;

PROC TTEST DATA=algae1; /* Testing differences between means using T-test */
CLASS group;
VAR chem3;
TITLE 'T test for significant difference in mean chem3 value for rivers of small/medium and large sizes';
RUN;


